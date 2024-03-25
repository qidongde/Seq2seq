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
batch_size = 16
hidden_size = 50
epochs = 20
mylr = 1e-4
dropout_p = 0.1


def train_test_split_func():
    train_pairs, test_pairs = []
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
        self.gru = nn.GRU(input_size, hidden_size, batch_first=True)

    def forward(self, input, hidden):
        # input[batch_size,time_step,input_size]
        # --> input_[batch_size,time_step,hidden_size]
        input_ = F.relu(self.embedding(input))

        # gru([batch_size,time_step,hidden_size],[num_layers,batch_size,hidden_size])
        # --> [batch_size,time_step,hidden_size] [num_layers,batch_size,hidden_size]
        output, hidden = self.gru(input_, hidden)
        return output, hidden

    def inithidden(self):
        # h0: [num_layers,batch_size,hidden_size]
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


class AttnDecoderRNN(nn.Module):
    def __init__(self, input_size, hidden_size, seq_length, dropout_p=0.1):
        super(AttnDecoderRNN, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.seq_length = seq_length
        self.dropout_p = dropout_p

        self.embedding = nn.Linear(self.input_size, self.hidden_size)
        self.attn = nn.Linear(self.hidden_size * 2, self.seq_length)
        self.attn1_combine = nn.Linear(self.hidden_size * 2, self.hidden_size)
        self.dropout = nn.Dropout(self.dropout_p)
        self.gru = nn.GRU(self.hidden_size, self.hidden_size, batch_first=True)

        self.attn2 = nn.Linear(self.hidden_size * 2, self.hidden_size)
        self.softmax = nn.LogSoftmax(dim=-1)

    def forward(self, input, hidden, encoder_outputs):
        # input(q): [batch_size,1,input_size]
        # hidden(k): [1,batch_size,hidden_size]
        # encoder_outputs(k): [batch_size,seq_length,hidden_size]

        # embedded: [batch_size,1,input_size] --> [batch_size,1,hidden_size]
        embedded = F.relu(self.embedding(input))
        embedded = self.dropout(embedded)

        # attn1_weights: [batch_size,seq_length]
        attn1_weights = F.softmax(
            self.attn(torch.cat((embedded.squeeze(), hidden.squeeze()), 1)), dim=1)

        # attn_applied[batch_size,1,hidden_size]
        # [batch_size,1,seq_length],[batch_size,seq_length,hidden_size] ---> [batch_size,1,hidden_size]
        attn_applied = torch.bmm(attn1_weights.unsqueeze(1), encoder_outputs)

        # weighted_input[batch_size,1,hidden_size]
        weighted_input = torch.cat((embedded.squeeze(), attn_applied.squeeze()), 1)
        weighted_input = self.attn_combine(weighted_input).unsqueeze(0)
        weighted_input = F.relu(weighted_input)

        # output: [batch_size,1,hidden_size]
        # hidden: [1,batch_size,hidden_size]
        output, hidden = self.gru(weighted_input, hidden)

        # attn2_weights: [batch_size,hidden_size]
        attn2_weights = F.softmax(
            self.attn2(torch.cat((embedded.squeeze(), encoder_outputs[:, -1, :]), 1)), dim=1)

        # [batch_size,1,hidden_size],[batch_size,hidden_size,1] ---> [batch_size,1,1]
        output = torch.bmm(output, attn2_weights.unsqueeze(-1)).squeeze()

        return output, hidden, attn1_weights, attn2_weights

    def inithidden(self):
        return torch.zeros(1, batch_size, self.hidden_size, device=device)


def Train_Iters(x, y, my_encoderrnn, my_attndecoderrnn, myadam_encode, myadam_decode, mse_loss):
    # encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
    encode_hidden = my_encoderrnn.inithidden()
    encode_output, encode_hidden = my_encoderrnn(x, encode_hidden)
    # [batch_size,121,25],[1,batch_size,hidden_size]
    # --> [batch_size,121,hidden_size],[1,batch_size,hidden_size]
    # encode_output_c [batch_size,121,hidden_size]

    decode_hidden = encode_hidden

    myloss = 0.0
    input_y = torch.zeros((batch_size, 1, decode_hidden), device=device)
    y_len = y.shape[1]
    y_pre = torch.zeros(batch_size, y_len, device=device)
    for idx in range(y_len):
        output_y, decode_hidden, attn_weight = my_attndecoderrnn(input_y, decode_hidden, encode_output)
        y_pre[:, idx] = output_y
        y_true = y[:, idx]
        myloss += mse_loss(output_y, y_true)
        input_y = output_y.detach()

    myadam_encode.zero_grad()
    myadam_decode.zero_grad()

    myloss.backward()

    myadam_encode.step()
    myadam_decode.step()

    return y_pre, y


def Train_seq2seq():
    train_pairs, test_pairs = train_test_split_func()

    train_dataset = MyPairsDataset(train_pairs)
    test_dataset = MyPairsDataset(test_pairs)
    train_dataloader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    test_dataloader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=False, drop_last=False)

    my_encoderrnn = EncoderRNN(25, hidden_size)
    my_attndecoderrnn = AttnDecoderRNN(hidden_size, hidden_size, 121, dropout_p=dropout_p)
    myadam_encode = optim.Adam(my_encoderrnn.parameters(), lr=mylr)
    myadam_decode = optim.Adam(my_attndecoderrnn.parameters(), lr=mylr)

    mse_loss = nn.MSELoss(reduction='sum')

    for epoch_idx in range(1, epochs + 1):
        train_y_true = []
        train_y_pre = []
        test_y_true = []
        test_y_pre = []

        my_encoderrnn.train()
        my_attndecoderrnn.train()
        for train_item, (train_x, train_y) in enumerate(tqdm(train_dataloader), start=1):
            y_pre, y = Train_Iters(train_x, train_y, my_encoderrnn, my_attndecoderrnn, myadam_encode, myadam_decode, mse_loss)

            train_y_true.extend(y.squeeze().tolist())
            train_y_pre.extend(y_pre.squeeze().tolist())

        eval_loss = nn.MSELoss()
        train_rmse_loss = np.sqrt(eval_loss(train_y_true, train_y_pre))
        
        print(f'The result of epoch{epoch_idx}:')
        print("Train RMSELoss:", train_rmse_loss)
        print("*" * 50)


if __name__ == '__main__':
    # dm_test_MyPairsDataset()
    Train_seq2seq()
