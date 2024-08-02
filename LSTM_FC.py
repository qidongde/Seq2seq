import mat73
from scipy.io import savemat
import numpy as np

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

from sklearn.metrics import mean_squared_error
import math
import matplotlib.pyplot as plt
import time

start_time = time.time()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device("cpu")
batch_size = 8
hidden_size = 200
dropout_ratio = 0.4
num_layers = 1
mylr = 1e-4
epochs = 20


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


class LSTM_FC(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, num_layers, batch_size):
        super(LSTM_FC, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.batch_size = batch_size

        self.lstm = nn.LSTM(input_size=input_size, hidden_size=self.hidden_size, num_layers=num_layers,
                            batch_first=True)
        self.linear1 = nn.Linear(self.hidden_size, 99)
        self.linear2 = nn.Linear(99, 99)
        self.tanh = nn.Tanh()
        self.dropout = nn.Dropout(p=dropout_ratio)

    def forward(self, input, hidden, c):
        # input(batch_size, seq_length, input_size)
        rr, (hn, c) = self.lstm(input, (hidden, c))
        # rr(batch_size, seq_len, hidden_size)
        # hn(num_layers, batch_size, hidden_size)
        # cn(num_layers, batch_size, hidden_size)
        output = rr[:, -1, :]
        output = self.dropout(output)
        output = self.linear1(output)
        output = self.linear2(self.tanh(output)).squeeze()

        output = output.to(device)
        hn = hn.to(device)
        c = c.to(device)

        return output, hn, c

    def inithiddenAndC(self):
        c0 = torch.zeros(num_layers, self.batch_size, self.hidden_size, device=device)
        h0 = torch.zeros(num_layers, self.batch_size, self.hidden_size, device=device)
        return h0, c0


def LSTM_train():
    train_pairs, test_pairs = train_test_split_func()

    train_dataset = MyPairsDataset(train_pairs)
    test_dataset = MyPairsDataset(test_pairs)
    train_dataloader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    test_dataloader = DataLoader(dataset=test_dataset, batch_size=batch_size, shuffle=True, drop_last=True)

    LSTM_model = LSTM_FC(22, hidden_size, 99, num_layers, batch_size).to(device)
    myadam = torch.optim.Adam(LSTM_model.parameters(), lr=mylr)
    mse_loss = nn.MSELoss()

    train_rmse_loss_final = 1e10
    test_rmse_loss_final = 1e10
    train_loss_list = []
    test_loss_list = []

    for epoch_idx in range(1, epochs + 1):
        train_y_true = []
        train_y_pre = []
        test_y_true = []
        test_y_pre = []
        mdic = {}

        LSTM_model.train()
        for train_item, (train_x, train_y) in enumerate(tqdm(train_dataloader), start=1):
            train_h0, train_c0 = LSTM_model.inithiddenAndC()
            train_output, train_hidden, train_c = LSTM_model(train_x, train_h0, train_c0)
            # print(type(train_output))
            # print(train_output.size())

            myadam.zero_grad()
            train_loss = mse_loss(train_output, train_y)
            train_loss.backward()
            myadam.step()

            train_y_true += train_y.tolist()
            train_y_pre += train_output.tolist()

        LSTM_model.eval()
        with torch.no_grad():
            for test_item, (test_x, test_y) in enumerate(test_dataloader, start=1):
                # print('test_x', test_x[0])
                # print('test_y', test_y[0])
                test_h0, test_c0 = LSTM_model.inithiddenAndC()
                test_predict, test_hidden, test_c = LSTM_model(test_x, test_h0, test_c0)

                test_y_true += test_y.tolist()
                test_y_pre += test_predict.tolist()
                # print('test_y_true:', test_y_true[-1])
                # print('test_y_pre:', test_y_pre[-1])

        train_rmse_loss = np.sqrt(mean_squared_error(train_y_true, train_y_pre))
        train_loss_list.append(train_rmse_loss)
        test_rmse_loss = np.sqrt(mean_squared_error(test_y_true, test_y_pre))
        test_loss_list.append(test_rmse_loss)

        print(f'The result of epoch{epoch_idx}:')
        print("Test RMSELoss:", test_rmse_loss)
        # print('Test R2', test_r2score)
        print("Train RMSELoss:", train_rmse_loss)
        # print("Train R2:", train_r2score)
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

            # np.save('output_save/LSTM_test_y_pre.npy', np.array(test_y_pre_final))
            torch.save(LSTM_model.state_dict(), './LSTM_FC_model_save/LSTM_FC_200_%d.pth' % epoch_idx)
            savemat("LSTM_FC_200_Result.mat", mdic)

    plt.figure()
    plt.plot(list(range(1, epochs + 1)), train_loss_list, label='train_loss', color='b')
    plt.plot(list(range(1, epochs + 1)), test_loss_list, label='test_loss', color='g')
    plt.legend(loc='best')
    plt.xlabel('Epoch')
    plt.ylabel('RMSE Loss')
    plt.xticks(np.arange(1, epochs + 1, 1))
    plt.savefig('./LSTM_FC_200_RMSE_loss.png')
    plt.show()

    end_time = time.time()
    time_consuming = end_time - start_time
    print(f'train_rmse_loss_final:{train_rmse_loss_final:.2f}')
    print(f'test_rmse_loss_final:{test_rmse_loss_final:.2f}')
    print(f'Time consuming:{time_consuming:.2f}s')


if __name__ == '__main__':
    LSTM_train()
