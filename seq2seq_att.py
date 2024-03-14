import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import torch.optim as optim
import time

import random
import matplotlib.pyplot as plt

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
batch_size = 16


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

        self.gru = nn.GRU(input_size, hidden_size, batch_first=True)

    def forward(self, input, hidden):
        # input[batch_size,time_step,input_size]
        output = self.embedding(input)

        # gru([batch_size,time_step,input_size],[num_layers,batch_size,hidden_size])
        # --> [batch_size,time_step,input_size] [num_layers,batch_size,hidden_size]
        output, hidden = self.gru(output, hidden)
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
        embedded = self.embedding(input)
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


if __name__ == '__main__':
    dm_test_MyPairsDataset()
