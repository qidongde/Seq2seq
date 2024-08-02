import math

import mat73
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.io import savemat
from sklearn.metrics import mean_squared_error
from torch.utils.data import Dataset, DataLoader

import torch.optim as optim
from tqdm import tqdm
import time

import random
import matplotlib.pyplot as plt

start_time = time.time()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
batch_size = 8
hidden_size = 50

output_size = 1
epochs = 20
mylr = 5e-5
dropout_p = 0.4
teacher_forcing_ratio = 0.5


def train_test_split_func():
    raw_data = mat73.loadmat('./Dataset_input_output_2014_2023.mat')
    # time_string = raw_data['time_selected']
    time_vector = raw_data['time_selected_vector']
    input_ = np.array(raw_data['input_selected']).squeeze()
    sigma = np.std(input_, axis=0)
    sigma[sigma == 0] = 1
    input = (input_ - np.mean(input_, axis=0)) / sigma
    output = raw_data['output_selected']

    test_filter = (time_vector[:, 0] == 2022)
    train_filter = (time_vector[:, 0] != 2022)

    x_train = np.transpose(input[train_filter, :, :], (0, 2, 1))
    y_train = np.log10(output[train_filter] + 1)
    x_test = np.transpose(input[test_filter, :, :], (0, 2, 1))
    y_test = np.log10(output[test_filter] + 1)

    train_pairs = list(zip(x_train, y_train))
    test_pairs = list(zip(x_test, y_test))
    # print('y_train.shape', y_train.shape)
    # print('y_test.shape', y_test.shape)
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
        # index = min(max(index, 0), self.sample_len - 1)

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
        input_ = torch.tanh(self.embedding(input))

        # gru([batch_size,time_step,hidden_size],[num_layers,batch_size,hidden_size])
        # --> [batch_size,time_step,hidden_size] [num_layers,batch_size,hidden_size]
        output, hidden = self.gru(input_, hidden)
        return output.to(device), hidden.to(device)

    def inithidden(self):
        # h0: [num_layers,batch_size,hidden_size]
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


class AttnDecoderRNN(nn.Module):
    def __init__(self, output_size, hidden_size, dropout_p=0.1):
        # output_size: 1
        # hidden_size: 200
        # dropout_p: 0.1
        # max_length: 99
        super(AttnDecoderRNN, self).__init__()
        self.output_size = output_size
        self.hidden_size = hidden_size
        self.dropout_p = dropout_p

        self.embedding = nn.Linear(1, hidden_size)
        self.attn = nn.Linear(self.hidden_size * 2, 121)
        self.attn_combine = nn.Linear(self.hidden_size * 2, self.hidden_size)

        self.dropout = nn.Dropout(self.dropout_p)
        self.gru = nn.GRU(self.hidden_size, self.hidden_size, batch_first=True)

        # (200,1)
        self.out = nn.Linear(self.hidden_size, self.output_size)

    def forward(self, input, hidden, encoder_outputs):
        # input(q): [batch_size,1]
        # hidden(k): [1,batch_size,200]
        # encoder_outputs(v): [batch_size,121,200]

        # [batch_size,1] --> [batch_size,200]
        embedded = self.embedding(input)

        # avoid overfitting
        embedded = self.dropout(embedded)

        # 1 attn_weights[batch_size,121]
        attn_weights = F.softmax(
            self.attn(torch.cat((embedded, hidden[0, :, :]), 1)), dim=-1)

        # 2 attn_applied[batch_size,1,200]
        # [batch_size,1,121],[batch_size,121,200] ---> [batch_size,1,200]
        attn_applied = torch.bmm(attn_weights.unsqueeze(1), encoder_outputs)

        # 3 output[batch_size,1,200]
        output = torch.cat((embedded, attn_applied[:, 0, :]), 1)
        output = self.attn_combine(output).unsqueeze(1)

        output = torch.tanh(output)

        # [batch_size,1,200],[1,batch_size,200] --> [batch_size,1,200],[1,batch_size,200]
        output, hidden = self.gru(output, hidden)
        # [batch_size,1,200]->[batch_size,200]->[batch_size,1]
        output = self.out(output[:, 0, :])

        # output[batch_size,1] hidden[1,batch_size,256] attn_weights[1,121]
        return output.to(device), hidden.to(device), attn_weights.to(device)

    def inithidden(self):
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


def Train_Iters(x, y, my_encoderrnn, my_attndecoderrnn, myadam_encode, myadam_decode, mymseloss):
    # 1 encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
    encode_hidden = my_encoderrnn.inithidden()
    # [batch_size,121,22],[1,batch_size,200] --> [batch_size,121,200] [1,batch_size,200]
    encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)

    # encode_output_c [batch_size,121,200]
    encode_output_c = encode_output
    # decode_hidden [1,batch_size,200]
    decode_hidden = encode_hidden

    input_y = torch.zeros((batch_size, 1), device=device)
    y_pre = torch.zeros((batch_size, 99), device=device)

    myloss = 0.0
    y_len = y.shape[1]

    use_teacher_forcing = True if random.random() < teacher_forcing_ratio else False
    if use_teacher_forcing:
        for idx in range(y_len):
            # [batch_size,1],[1,batch_size,200],[batch_size,121,200] ---> [batch_size,1],[1,batch_size,200],[1,121]
            output_y, decode_hidden, attn_weight = my_attndecoderrnn(input_y, decode_hidden, encode_output_c)
            y_pre[:, idx] = output_y.squeeze()
            target_y = y[:, idx].unsqueeze(1)
            myloss = myloss + mymseloss(output_y, target_y)
            input_y = target_y
    else:
        for idx in range(y_len):
            # [batch_size,1],[1,batch_size,200],[batch_size,121,200] ---> [batch_size,1],[1,batch_size,200],[1,121]
            output_y, decode_hidden, attn_weight = my_attndecoderrnn(input_y, decode_hidden, encode_output_c)
            y_pre[:, idx] = output_y.squeeze()
            target_y = y[:, idx].unsqueeze(1)
            myloss = myloss + mymseloss(output_y, target_y)
            input_y = output_y.detach()

    myadam_encode.zero_grad()
    myadam_decode.zero_grad()

    myloss.backward()

    myadam_encode.step()
    myadam_decode.step()

    return y, y_pre, myloss.item() / y_len


def Seq2Seq_Evaluate(x, my_encoderrnn, my_attndecoderrnn):
    with torch.no_grad():
        # 1 encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
        encode_hidden = my_encoderrnn.inithidden()
        # [1,121,22],[1,1,200] --> [1,121,200] [1,1,200]
        encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)

        # encode_output_c [1,121,200]
        encoder_outputs_c = encode_output
        # decode_hidden [1,1,200]
        decode_hidden = encode_hidden

        input_y = torch.zeros((batch_size, 1), device=device)
        y_pre = torch.zeros((batch_size, 99), device=device)
        decoder_attentions = torch.zeros((batch_size, 99, 121), device=device)

        for idx in range(99):
            output_y, decode_hidden, attn_weights = my_attndecoderrnn(input_y, decode_hidden, encoder_outputs_c)
            y_pre[:, idx] = output_y.squeeze()
            decoder_attentions[:, idx, :] = attn_weights
            input_y = output_y.detach()

        return y_pre, decoder_attentions


def Train_seq2seq():
    train_pairs, test_pairs = train_test_split_func()
    train_dataset = MyPairsDataset(train_pairs)
    test_dataset = MyPairsDataset(test_pairs)
    train_dataloader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    test_dataloader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    my_encoderrnn = EncoderRNN(22, hidden_size).to(device)
    my_attndecoderrnn = AttnDecoderRNN(output_size=output_size, hidden_size=hidden_size, dropout_p=dropout_p).to(device)

    myadam_encode = optim.Adam(my_encoderrnn.parameters(), lr=mylr)
    myadam_decode = optim.Adam(my_attndecoderrnn.parameters(), lr=mylr)

    mymseloss = nn.MSELoss()

    train_rmse_loss_final = 1e10
    test_rmse_loss_final = 1e10
    train_loss_list = []
    test_loss_list = []

    for epoch_idx in range(1, 1 + epochs):

        train_y_true = []
        train_y_pre = []
        test_y_true = []
        test_y_pre = []
        mdic = {}

        my_encoderrnn.train()
        my_attndecoderrnn.train()
        for item, (x, y) in enumerate(tqdm(train_dataloader), start=1):
            train_y, train_output, myloss = Train_Iters(x, y, my_encoderrnn, my_attndecoderrnn, myadam_encode,
                                                        myadam_decode,
                                                        mymseloss)
            train_y_true += train_y.tolist()
            train_y_pre += train_output.tolist()

        my_encoderrnn.eval()
        my_attndecoderrnn.eval()
        with torch.no_grad():
            for test_item, (test_x, test_y) in enumerate(test_dataloader, start=1):
                test_predict, decoder_attentions = Seq2Seq_Evaluate(test_x, my_encoderrnn, my_attndecoderrnn)
                test_y_true += test_y.tolist()
                test_y_pre += test_predict.tolist()

        train_rmse_loss = np.sqrt(mean_squared_error(train_y_true, train_y_pre))
        train_loss_list.append(train_rmse_loss)
        test_rmse_loss = np.sqrt(mean_squared_error(test_y_true, test_y_pre))
        test_loss_list.append(test_rmse_loss)

        print(f'The result of epoch{epoch_idx}:')
        print("Test RMSELoss:", test_rmse_loss)
        print("Train RMSELoss:", train_rmse_loss)
        print("*" * 50)
        if test_rmse_loss < test_rmse_loss_final:
            train_rmse_loss_final = train_rmse_loss
            test_rmse_loss_final = test_rmse_loss
            train_y_true_final = train_y_true
            train_y_pre_final = train_y_pre
            test_y_true_final = test_y_true
            test_y_pre_final = test_y_pre
            mdic = {
                'train_y_true_final': train_y_true_final,
                'train_y_pre_final': train_y_pre_final,
                'test_y_true_final': test_y_true_final,
                'test_y_pre_final': test_y_pre_final
            }

            torch.save(my_encoderrnn.state_dict(), './seq2seq_model_save/my_encoderrnn_50_2_%d.pth' % epoch_idx)
            torch.save(my_attndecoderrnn.state_dict(), './seq2seq_model_save/my_attndecoderrnn_50_2_%d.pth' % epoch_idx)
            # savemat("Seq2Seq_50_Result.mat", mdic)

    plt.figure()
    plt.plot(list(range(1, epochs + 1)), train_loss_list, label='train_loss', color='b')
    plt.plot(list(range(1, epochs + 1)), test_loss_list, label='test_loss', color='g')
    plt.legend(loc='best')
    plt.xlabel('Epoch')
    plt.ylabel('RMSE Loss')
    plt.xticks(np.arange(1, epochs + 1, 1))
    plt.savefig('./seq2seq_50_2_RMSE_loss.png')
    plt.show()

    end_time = time.time()
    time_consuming = end_time - start_time
    print(f'train_rmse_loss_final:{train_rmse_loss_final:.2f}')
    print(f'test_rmse_loss_final:{test_rmse_loss_final:.2f}')
    print(f'Time consuming:{time_consuming:.2f}s')



def dm_test_Attention():
    train_pairs, test_pairs = train_test_split_func()
    train_dataset = MyPairsDataset(train_pairs)
    test_dataset = MyPairsDataset(test_pairs)
    train_dataloader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=False, drop_last=False)
    test_dataloader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=False, drop_last=False)

    my_encoderrnn = EncoderRNN(22, hidden_size).to(device)
    my_attndecoderrnn = AttnDecoderRNN(output_size=output_size, hidden_size=hidden_size, dropout_p=dropout_p).to(
        device)

    my_encoderrnn.load_state_dict(torch.load('seq2seq_model_save/my_encoderrnn_100_19.pth'))
    my_attndecoderrnn.load_state_dict(torch.load('seq2seq_model_save/my_attndecoderrnn_100_19.pth'))

    train_y_true = []
    train_y_pre = []
    test_y_true = []
    test_y_pre = []

    my_encoderrnn.eval()
    my_attndecoderrnn.eval()
    with torch.no_grad():
        for train_item, (train_x, train_y) in enumerate(tqdm(train_dataloader), start=1):
            train_predict, train_decoder_attentions = Seq2Seq_Evaluate(train_x, my_encoderrnn, my_attndecoderrnn)
            train_y_true += train_y.tolist()
            train_y_pre += train_predict.tolist()

        for test_item, (test_x, test_y) in enumerate(tqdm(test_dataloader), start=1):
            test_predict, test_decoder_attentions = Seq2Seq_Evaluate(test_x, my_encoderrnn, my_attndecoderrnn)
            test_y_true += test_y.tolist()
            test_y_pre += test_predict.tolist()

        mdic = {
            'train_y_true_final': train_y_true,
            'train_y_pre_final': train_y_pre,
            'test_y_true_final': test_y_true,
            'test_y_pre_final': test_y_pre
        }

        savemat("Seq2Seq_100_Result_.mat", mdic)


if __name__ == '__main__':
    # dm_test_MyPairsDataset()
    Train_seq2seq()
    # dm_test_Attention()
