import pandas as pd
import numpy as np


"""
                                                         Structure of the data

                              X1_(pos1-10) X1_(pos1-9) ... X2_(pos1-10) X2_(pos1-9) ... X3_...            y_(pos1)
                              X1_(pos2-10) X1_(pos2-9) ... X2_(pos2-10) X2_(pos2-9) ... X3_...            y_(pos2)
                              X1_(pos3-10) X1_(pos3-9) ... X2_(pos3-10) X2_(pos3-9) ... X3_...            y_(pos3)
                                                       ...

"""

"""
                          IMPORTANT !!!!
            The column order of the input data must match the order as returned by the make_data_optimized.py function

"""


def make_training_data(filtered_data, fname_X, fname_alt_counts, fname_ref_counts, windowsize_left = 10, windowsize_right = 10):
    
    data = pd.read_csv(filtered_data)
    X_list = []
    alt_counts_list = []
    ref_counts_list = []

    for i in range(windowsize_left+1, data.shape[0]-(windowsize_left+1)):

        print(str(i/data.shape[0]*100)+'%')
        
        alt_counts_list.append(int(data.loc[i,'general_alt_counts']))
        ref_counts_list.append(int(data.loc[i,'ref_counts']))
        
        # the predictors upstream and downstream of the given locus
        context_upstream = data.iloc[(i-windowsize_left):i,5:]
        context_downstream = data.iloc[(i+1):(i+windowsize_right),5:]
        x = np.concatenate([context_upstream.values.flatten(), context_downstream.values.flatten()])
        X_list.append(x)


    X = pd.DataFrame(X_list)
    alt_counts_list = pd.Series(alt_counts_list)
    ref_counts_list = pd.Series(ref_counts_list)

    X.to_csv(fname_X, index = False)
    alt_counts_list.to_csv(fname_alt_counts, index = False)
    ref_counts_list.to_csv(fname_ref_counts, index = False)
