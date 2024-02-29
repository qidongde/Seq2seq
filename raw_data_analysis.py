from scipy.io import loadmat
import mat73
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy.stats import norm

warnings.filterwarnings('ignore')

filename1 = 'Multi_year_CCN_LSTM_new_10days_newCCNdata_v3.mat'
input_data = mat73.loadmat(filename1)
all_features = np.array(input_data['data_input_valid']).squeeze()
all_features_flat = all_features.reshape(-1, 27)
feature_name_list = ['lat', 'lon', 'Height', 'Pressure', 'Land_flag', 'PBLH', 'PRECTOTCORR', 'CLDLOW', 'CLDMID',
                     'SWGDN', 'SLP', 'TQL', 'WS10M', 'LTS', 'TS', 'SO2', 'DU001', 'SS001', 'OMEGA', 'QL', 'RH', 'T',
                     'EPV', 'CO', 'chla', 'DMSflux', 'DMSconc']


def feature_profile_analysis(num):
    feature_selected = all_features_flat[:, num]
    # feature_selected = np.where(feature_selected != 0, np.log(feature_selected), 0)
    plt.figure(figsize=(8, 4))
    sns.distplot(feature_selected, bins=100, hist=True, kde=False, norm_hist=False,
                 rug=False, vertical=False, label='Distplot',
                 axlabel=feature_name_list[num], hist_kws={'color': 'y', 'edgecolor': 'k'},
                 fit=norm)

    plt.legend()
    plt.grid(linestyle='--')
    plt.show()


if __name__ == '__main__':
    feature_profile_analysis(26)
