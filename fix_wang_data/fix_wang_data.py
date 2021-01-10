#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from AnalysisCode.global_variables import phenotypic_variables, create_folder, dataset_names, sm_datasets, mm_datasets, wang_datasets, tanouchi_datasets


checklist = pd.DataFrame(columns = ['dataset', 'filename', 'start', 'end', 'useful'])


for data_origin in wang_datasets:  # input_args.dataset_names:
    print(data_origin)
    
    temp_dict = {}
    temp_ind = 0
    
    """
    data_origin ==> Name of the dataset we are analysing
    raw_data ==> Where the folder containing the raw data for this dataset is
    processed_data ==> The folder we will put the processed data in
    """
    args = {
        'data_origin': data_origin,
        'raw_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/RawData/',
        'processed_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/'
    }

    infiles = glob.glob(args['raw_data'] + '/*.dat')

    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(infiles[5:]):
    
        filename = file.split('/')[-1].split('.')[0]
        
        # checklist = checklist.append({
        #     'dataset': data_origin,
        #     'filename': filename,
        #     'start': '',
        #     'end': '',
        #     'useful': ''
        # }, ignore_index=True)
        #
        # continue
        
        print(filename)

        raw_lineage = pd.read_csv(file, delimiter=' ')  # names=['time', 'division_flag', 'length', 'width', 'area', 'yfp_intensity', 'CMx', 'CMy']

        # Sometimes the time indices are called time and others called index
        if 'time' in raw_lineage.columns:
            raw_lineage = raw_lineage.rename(columns={'division': 'division_flag'})[['time', 'division_flag', 'length']].sort_values('time')
        elif 'index' in raw_lineage.columns:
            raw_lineage = raw_lineage.rename(columns={'index': 'time', 'division': 'division_flag'})[['time', 'division_flag', 'length']].sort_values('time')

        raw_lineage['length'] = (raw_lineage['length']) * 0.0645  # Convert it from pixel length to micrometers
        
        print('end:', len(raw_lineage)-1)
        
        fig, ax = plt.subplots(figsize=[13, 7])
        
        # Show the series
        plt.plot(raw_lineage['time'].values, raw_lineage['length'].values)
        plt.title(filename)
        plt.xlabel('index')
        plt.ylabel('length (micrometers)')
        plt.tight_layout()
        plt.show()
        plt.close()
        
        continue

        useful = np.int(input('useful'))
        if useful:
            start = np.int(input('start'))
            print('start is {}'.format(start))
            end = np.int(input('end'))
            print('end is {}'.format(end))
        else:
            start = np.nan
            end = np.nan
        
        if useful:
            plt.plot(raw_lineage['time'].iloc[start:end].values, raw_lineage['length'].iloc[start:end].values)
            plt.xlabel('index')
            plt.ylabel('length (micrometers)')
            plt.tight_layout()
            plt.show()
            plt.close()
        
        temp_dict[count] = {
            'dataset': data_origin,
            'filename': filename,
            'start': start,
            'end': end,
            'useful': useful
        }

    # print('temp_dict:', temp_dict, sep='\n')
    # # Append them to the big dataframes
    # checklist = checklist.append(pd.DataFrame.from_dict(temp_dict, "index"), ignore_index=True)
    # print('checklist:', checklist, sep='\n')
    
    print('finished {}!!!'.format(data_origin))

# save it
checklist.to_csv('checklist.csv', index=False)


