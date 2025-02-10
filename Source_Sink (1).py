#%%
# Import Packages
##########################
import numpy as np
import seaborn as sns
from ieeg.auth import Session
from sklearn.linear_model import LinearRegression
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import os
import re
import sys
import seaborn as sns
from scipy.signal import butter, filtfilt
import scipy.signal as sc
from utils import get_iEEG_data, plot_iEEG_data, artifact_removal, notch_filter, bandpass_filter, stopband_filter, detect_bad_channels # remove for scalp!


######### Functions ######### 
# Preprocessing (change for scalp!)
def preprocess_data(dataset_name,start_time, pw_path,usr_name):

    num_mins = 2 # Use 2 or 5 mins usually
    full_start_time = start_time
    full_stop_time = full_start_time + float(num_mins*60)

    full_start_usec = full_start_time * 1e6
    full_stop_usec = full_stop_time * 1e6

    df, fs = get_iEEG_data(
        username=usr_name,
        password_bin_file=pw_path,
        iEEG_filename=dataset_name,
        start_time_usec=full_start_usec,
        stop_time_usec=full_stop_usec
    )

    channel_names = df.columns.tolist() 
    data = df.to_numpy()

    # Apply bandpass filter
    data = bandpass_filter(data, fs, order=4, lo=0.5, hi=150)

    # Apply notch filters at 60 Hz and harmonics
    data = notch_filter(data, fs)

    # Detect bad channels for iEEG 
    channel_mask, details = detect_bad_channels(data, fs)

    # Filter out bad channels
    data = data[:, channel_mask]
    valid_channels = [ch for i, ch in enumerate(channel_names) if channel_mask[i]]

    # Common Average Reference
    data -= np.mean(data, axis=0)
    
    # Downsample 
    target_fs = 256
    fs = int(fs)
    # Calculate the resampling factors
    gcd = np.gcd(fs, target_fs)
    up = target_fs // gcd
    down = fs // gcd

    # Resample each channel using polynomial interpolation
    data = sc.resample_poly(data, up, down, axis=0)

    # New sampling frequency
    fsd = target_fs

    return data, valid_channels, fsd

# Fit LTI Model on 500 ms Interictal Data : x(t + 1) = Ax(t)
def fit_ar_model(data, fs, window_length=0.5):

    # 500 ms windows chosen in paper
    window_size = window_length*fs
    n_windows = data.shape[0] // window_size
    model = LinearRegression()
    
    A_all = []
    for i in range(int(n_windows)):
        window_data = data[int(i * window_size):int((i + 1) * window_size), :]
        X = window_data[:-1]
        y = window_data[1:]
        model.fit(X, y)
        A = model.coef_
        A_all.append(A)
    
    A_all = np.array(A_all)
    Avg_A = np.mean(A_all, axis=0)
    Avg_A = np.abs(Avg_A)
    print(f'Avg A: {Avg_A}')
    print(f'Avg A shape : {Avg_A.shape}')

    # Avg A is a matrix of size NxN , where N is number of channels!
    # i = influenced
    # j = influences!

    return(Avg_A)

def calculate_rank(arr): 
    # Get the indices that would sort the array in descending order
    sorted_indices = np.argsort(arr)[::-1]
    # Create a ranking array
    ranks = np.zeros_like(arr, dtype=float)
    ranks[sorted_indices] = np.arange(1, len(arr) + 1) / len(ranks)
    print(f'rank = {ranks}')
    return ranks

# Source Sink Space
def calculate_source_sink_space(Avg_A):
    np.fill_diagonal(Avg_A, 0)# Ignore self
    i_node_strength = np.sum(Avg_A, axis=0)
    j_node_strength = np.sum(Avg_A, axis=1)

    # rank node strengths to get row rank (rr) and column rank (cr)
    rr = calculate_rank(i_node_strength)
    cr = calculate_rank(j_node_strength)
    print(f'rr = {rr}')
    print(f'cr = {cr}')

    plt.scatter(rr, cr)
    plt.title('Source Sink Space')
    plt.xlabel('Influence by others (rr)')
    plt.ylabel('Influence to others (cr)')
    plt.close()
    return(rr, cr, i_node_strength, j_node_strength)

# Calculate Source Index
def calculate_source_index(rr,cr):
    N = len(rr)
    x = rr - (1/N)
    y = cr - 1
    vector_length = np.sqrt(x**2 + y**2)
    print(f'vector_length: {vector_length}')
    source_index = np.sqrt(2) - vector_length
    return(source_index)

# Calculate Sink Index
def calculate_sink_index(rr,cr):
    N = len(rr)
    x = rr - 1
    y = cr - (1/N)
    vector_length = np.sqrt(x**2 + y**2)
    sink_index = np.sqrt(2) - vector_length
    return(sink_index)

def calculate_source_influence_sink_connectivity(Avg_A, source_index, sink_index):

    # influence i = sum(A(ij)*source_index(j))
    source_influence = np.dot(Avg_A, source_index)
    sink_connectivity = np.dot(Avg_A, sink_index)

    return source_influence, sink_connectivity

# Calculate All
def calculate_ss_measures(filename, start_time, pw_path,usr_name):
    data, valid_channels, fsd = preprocess_data(filename,start_time, pw_path,usr_name)
    Avg_A = fit_ar_model(data, fsd, window_length=0.5)
    rr,cr,i_node_strength, j_node_strength = calculate_source_sink_space(Avg_A)
    source_index = calculate_source_index(rr,cr)
    sink_index = calculate_sink_index(rr,cr)

    source_influence, sink_connectivity = calculate_source_influence_sink_connectivity(Avg_A, source_index, sink_index)
    
    df = pd.DataFrame({
        'channel' : valid_channels,
        'source_index' : source_index,
        'sink_index' : sink_index,
        'source_influence' : source_influence,
        'sink connectivity': sink_connectivity
    })

    return(df)

#%%
# Input table with patient ID, filenames and starttimes 
table = pd.read_csv('your_patient_table.csv')

#iEEG recon password and username - For Example:
pw_path = '../../sla_ieeglogin.bin'
usr_name = 'slavelle'

#output directory
outputdir = ''

# filename, starttime
for _ ,row in table.iterrows():
    try:
        pt_id = row['HUPID']
        # filename = row['filename']
        # start_time = row['start_time']

        filename, start_time  = row['CCEP_baseline_training_5mins'].split(',')
        start_time = float(start_time)

        ss_df = calculate_ss_measures(filename, start_time, pw_path, usr_name)
        print(ss_df)

        ss_df.to_pickle(f'{outputdir}/{pt_id}_source_sink.pkl')
    except Exception as e:
        print(f'{pt_id} skipped: {e}')
