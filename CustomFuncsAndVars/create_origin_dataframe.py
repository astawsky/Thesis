#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, mm_data_names, symbols, seaborn_preamble, shuffle_info, create_folder
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from scipy.stats import linregress


origin = pd.DataFrame(columns=['file', 'folder(s)', '#_times_repeated', 'data_extension'])

folder_names = glob.glob('/Users/alestawsky/PycharmProjects/Thesis/RawData/MM/*')

for folder_full in folder_names:
    print(folder_full)
    file_names = glob.glob(folder_full + '/*')
    folder = folder_full.split('/')[-1]
    for file_full in file_names:
        
        file = file_full.split('/')[-1].split('.')[0]
        
        if file in origin.file.values:
            print('INFILE!')
            folders = origin[origin['file'] == file]['folder(s)'].values[0] + ', ' + folder
            # print(folders.split(', '))
            # print(len(folders.split(', ')))
            # print(type(folders))
            # exit()
            origin.loc[origin['file'] == file, 'folder(s)'] = folders
            origin.loc[origin['file'] == file, '#_times_repeated'] = len(folders.split(', '))
            origin.loc[origin['file'] == file, 'data_extension'] = file_full.split('/')[-1].split('.')[-1]
        else:
            
            origin = origin.append({
                'file': file,
                'folder(s)': str(folder),
                '#_times_repeated': 1,
                'data_extension': file_full.split('/')[-1].split('.')[-1]
            }, ignore_index=True)


# for col in origin.columns:
#     print(origin[col])

print(origin)

origin.to_csv('/Users/alestawsky/PycharmProjects/Thesis/RawData/repetitions_of_files.csv', index=False)

print(origin['#_times_repeated'])
print('len    ', len(origin['#_times_repeated'].unique()))
print(origin['#_times_repeated'].unique())
