import math

import mat73
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import torch.optim as optim
from tqdm import tqdm
import time

import random
import matplotlib.pyplot as plt

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
batch_size = 1
hidden_size = 256
# max dn/dlogp = 2590
output_size = 2950
epochs = 50
mylr = 1e-3
dropout_p = 0.1
print_interval_num = 50
plot_interval_num = 10
teacher_forcing_ratio = 0.5


def train_test_split_func():
    raw_data = mat73.loadmat('./Dataset_selected_2016.mat')
    time_string = raw_data['time_selected']
    # time_vector = raw_data['time_selected_vector']
    input = np.array(raw_data['input_selected']).squeeze()
    output = np.around(raw_data['output_selected'])

    split_label = np.floor(time_string) % 7
    test_filter = split_label == 3
    train_filter = split_label != 3

    x_train = np.transpose(input[train_filter, :, :], (0, 2, 1))
    y_train = output[train_filter]
    x_test = np.transpose(input[test_filter, :, :], (0, 2, 1))
    y_test = output[test_filter]

    train_pairs = list(zip(x_train, y_train))
    test_pairs = list(zip(x_test, y_test))

    return train_pairs, test_pairs


# raw_data -> MyPairsDataset --> DataLoader
# 1 __init__(self, my_pairs) self.my_pairs self.sample_len
# 2 __len__(self)
# 3 __getitem__(self, index)
# Dataset:[[x1,y1],[x2,y2]...]
class MyPairsDataset(Dataset):
    def __init__(self, my_pairs):
        self.my_pairs = my_pairs
        self.sample_len = len(my_pairs)

    def __len__(self):
        return self.sample_len

    def __getitem__(self, index):
        index = min(max(index, 0), self.sample_len - 1)

        x = self.my_pairs[index][0]
        y = self.my_pairs[index][1]

        tensor_x = torch.tensor(x, dtype=torch.float, device=device)
        tensor_y = torch.tensor(y, dtype=torch.float, device=device)
        # print('tensor_y.shape===>', tensor_y.shape, tensor_y)

        return tensor_x, tensor_y


def dm_test_MyPairsDataset():
    x_train = np.random.random((2000, 25, 121))
    y_train = np.random.random((2000, 30))
    my_pairs = list(zip(x_train, y_train))
    mypairsdataset = MyPairsDataset(my_pairs)
    mydataloader = DataLoader(dataset=mypairsdataset, batch_size=batch_size, shuffle=True)
    for i, (x, y) in enumerate(mydataloader):
        print('x.shape', x.shape)
        print('y.shape', y.shape)
        if i == 1:
            break


class EncoderRNN(nn.Module):
    def __init__(self, input_size, hidden_size):
        # input_size: feature number
        # hidden_size: hyperparameter
        super(EncoderRNN, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size

        self.embedding = nn.Linear(input_size, hidden_size)
        self.gru = nn.GRU(hidden_size, hidden_size, batch_first=True)

    def forward(self, input, hidden):
        # input[batch_size,time_step,input_size]
        # --> input_[batch_size,time_step,hidden_size]
        input_ = F.relu(self.embedding(input))

        # gru([batch_size,time_step,hidden_size],[num_layers,batch_size,hidden_size])
        # --> [batch_size,time_step,hidden_size] [num_layers,batch_size,hidden_size]
        output, hidden = self.gru(input_, hidden)
        return output.to(device), hidden.to(device)

    def inithidden(self):
        # h0: [num_layers,batch_size,hidden_size]
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


class AttnDecoderRNN(nn.Module):
    def __init__(self, output_size, hidden_size, dropout_p=0.1):
        # output_size: 2590
        # hidden_size: 256
        # dropout_p: 0.1
        # max_length: 99
        super(AttnDecoderRNN, self).__init__()
        self.output_size = output_size
        self.hidden_size = hidden_size
        self.dropout_p = dropout_p

        # nn.Embedding(2590,256)
        self.embedding = nn.Embedding(self.output_size, self.hidden_size)
        self.attn = nn.Linear(self.hidden_size * 2, 121)
        self.attn_combine = nn.Linear(self.hidden_size * 2, self.hidden_size)

        self.dropout = nn.Dropout(self.dropout_p)
        self.gru = nn.GRU(self.hidden_size, self.hidden_size, batch_first=True)

        # (256,2590)
        self.out = nn.Linear(self.hidden_size, self.output_size)

        # Normalization
        self.softmax = nn.LogSoftmax(dim=-1)

    def forward(self, input, hidden, encoder_outputs):
        # input(q): [batch_size,1]
        # hidden(k): [1,batch_size,256]
        # encoder_outputs(v): [batch_size,121,256]

        # [batch_size,1] --> [batch_size,1,256]
        embedded = self.embedding(input.long())

        # avoid overfitting
        embedded = self.dropout(embedded)

        # 1 attn_weights[batch_size,99]
        attn_weights = F.softmax(
            self.attn(torch.cat((embedded[:, 0, :], hidden[0, :, :]), 1)))

        # 2 attn_applied[1,1,256]
        # [batch_size,1,121],[batch_size,121,256] ---> [batch_size,1,256]
        attn_applied = torch.bmm(attn_weights.unsqueeze(1), encoder_outputs)

        # 3 output[batch_size,1,256]
        output = torch.cat((embedded[:, 0, :], attn_applied[:, 0, :]), 1)
        output = self.attn_combine(output).unsqueeze(1)

        output = F.relu(output)

        # [batch_size,1,256],[1,batch_size,256] --> [batch_size,1,256],[1,batch_size,256]
        output, hidden = self.gru(output, hidden)
        # [batch_size,1,256]->[batch_size,256]->[batch_size,2590]
        output = self.softmax(self.out(output[:, 0, :]))

        # output[batch_size,2590] hidden[1,batch_size,256] attn_weights[1,121]
        return output.to(device), hidden.to(device), attn_weights.to(device)

    def inithidden(self):
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


def Train_Iters(x, y, my_encoderrnn, my_attndecoderrnn, myadam_encode, myadam_decode, mycrossentropyloss):
    # 1 encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
    encode_hidden = my_encoderrnn.inithidden()
    # [batch_size,121,22],[1,batch_size,256] --> [batch_size,121,256] [1,batch_size,256]
    encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)

    # encode_output_c [batch_size,121,256]
    encode_output_c = encode_output
    # decode_hidden [1,batch_size,256]
    decode_hidden = encode_hidden

    input_y = torch.zeros((batch_size, 1), device=device)

    myloss = 0.0
    y_len = y.shape[1]

    use_teacher_forcing = True if random.random() < teacher_forcing_ratio else False
    if use_teacher_forcing:
        for idx in range(y_len):
            # [batch_size,1],[1,batch_size,256],[batch_size,121,256] ---> [batch_size,2950],[1,batch_size,256],[1,121]
            output_y, decode_hidden, attn_weight = my_attndecoderrnn(input_y, decode_hidden, encode_output_c)
            target_y = y[:, idx]
            myloss = myloss + mycrossentropyloss(output_y, target_y.long())
            input_y = target_y.unsqueeze(1)
    else:
        for idx in range(y_len):
            # [batch_size,1],[1,batch_size,256],[batch_size,121,256] ---> [batch_size,2950],[1,batch_size,256],[1,121]
            output_y, decode_hidden, attn_weight = my_attndecoderrnn(input_y, decode_hidden, encode_output_c)
            target_y = y[:, idx]
            myloss = myloss + mycrossentropyloss(output_y, target_y.long())
            topv, topi = output_y.topk(1)
            input_y = topi.detach()

    myadam_encode.zero_grad()
    myadam_decode.zero_grad()

    myloss.backward()

    myadam_encode.step()
    myadam_decode.step()

    return myloss.item() / y_len


def vec_cos_similarity(vec1, vec2):
    cos_similarity = np.sum(vec1 * vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))

    return cos_similarity


train_pairs, test_pairs = train_test_split_func()


def Train_seq2seq():
    train_dataset = MyPairsDataset(train_pairs)
    train_dataloader = DataLoader(dataset=train_dataset, batch_size=batch_size, drop_last=True, shuffle=True)

    my_encoderrnn = EncoderRNN(22, 256).to(device)
    my_attndecoderrnn = AttnDecoderRNN(output_size=2950, hidden_size=256, dropout_p=0.1).to(device)

    myadam_encode = optim.Adam(my_encoderrnn.parameters(), lr=mylr)
    myadam_decode = optim.Adam(my_attndecoderrnn.parameters(), lr=mylr)

    mycrossentropyloss = nn.NLLLoss()

    plot_loss_list = []

    for epoch_idx in range(1, 1 + epochs):

        print_loss_total, plot_loss_total = 0.0, 0.0
        starttime = time.time()

        for item, (x, y) in enumerate(tqdm(train_dataloader), start=1):
            myloss = Train_Iters(x, y, my_encoderrnn, my_attndecoderrnn, myadam_encode, myadam_decode,
                                 mycrossentropyloss)
            print_loss_total += myloss
            plot_loss_total += myloss

            if item % print_interval_num == 0:
                print_loss_avg = print_loss_total / print_interval_num
                print_loss_total = 0
                print('Epochs: %d  Loss: %.6f Time:%d' % (epoch_idx, print_loss_avg, time.time() - starttime))

            if item % plot_interval_num == 0:
                plot_loss_avg = plot_loss_total / plot_interval_num
                plot_loss_list.append(plot_loss_avg)
                plot_loss_total = 0

        torch.save(my_encoderrnn.state_dict(), './model_save/my_encoderrnn_%d.pth' % epoch_idx)
        torch.save(my_attndecoderrnn.state_dict(), './model_save/my_attndecoderrnn_%d.pth' % epoch_idx)

    plt.figure()
    plt.plot(plot_loss_list)
    plt.savefig('./s2sq_loss.png')
    plt.show()


def Seq2Seq_Evaluate(x, my_encoderrnn, my_attndecoderrnn):
    with torch.no_grad():
        # 1 encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
        encode_hidden = my_encoderrnn.inithidden()
        # [1,121,22],[1,1,256] --> [1,121,256] [1,1,256]
        encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)

        # encode_output_c [1,121,256]
        encoder_outputs_c = encode_output
        # decode_hidden [1,1,256]
        decode_hidden = encode_hidden

        input_y = torch.zeros((batch_size, 1), device=device)
        y_pre_list = []
        decoder_attentions = torch.zeros(99, 121)
        for idx in range(99):
            output_y, decode_hidden, attn_weights = my_attndecoderrnn(input_y, decode_hidden, encoder_outputs_c)
            topv, topi = output_y.topk(1)
            decoder_attentions[idx] = attn_weights
            input_y = topi.detach()
            y_pre_list.append(topi.item())
        return y_pre_list, decoder_attentions


PATH1 = './model_save/my_encoderrnn_48.pth'
PATH2 = './model_save/my_attndecoderrnn_48.pth'


def dm_test_Attention():
    testdataset = MyPairsDataset(test_pairs)
    testdataloader = DataLoader(dataset=testdataset, batch_size=1, shuffle=False)

    input_size = 22
    hidden_size = 256
    my_encoderrnn = EncoderRNN(input_size, hidden_size).to(device)
    my_encoderrnn.load_state_dict(torch.load(PATH1))
    # my_encoderrnn.load_state_dict(torch.load(PATH1, map_location=lambda storage, loc: storage), False)

    input_size = 2950
    hidden_size = 256
    my_attndecoderrnn = AttnDecoderRNN(input_size, hidden_size).to(device)
    my_attndecoderrnn.load_state_dict(torch.load(PATH2))
    # my_attndecoderrnn.load_state_dict(torch.load(PATH2, map_location=lambda storage, loc: storage), False)
    for item, (x, y) in enumerate(testdataloader, start=1):
        if item == 200:
            y_pre_list, decoder_attentions = Seq2Seq_Evaluate(x, my_encoderrnn, my_attndecoderrnn)
            print('y_pre -->', y_pre_list)
            print('y_true -->', y.tolist())

            plt.matshow(decoder_attentions.numpy())

            # plt.savefig("./s2s_attn.png")
            plt.show()
            # print('attentions.numpy()--->\n', attentions.numpy())
            # print('attentions.size--->', attentions.size())
            break


if __name__ == '__main__':
    # dm_test_MyPairsDataset()
    # Train_seq2seq()
    dm_test_Attention()
