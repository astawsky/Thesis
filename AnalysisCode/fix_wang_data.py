#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from AnalysisCode.global_variables import phenotypic_variables, create_folder, dataset_names, sm_datasets, mm_datasets, wang_datasets, tanouchi_datasets


def fix(args):
    
    # In case we can't use some files we want the lineage IDs to be in integer order
    offset = 0
    
    # The filenames
    infiles = glob.glob(args['raw_data'] + '/*.dat')

    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(infiles):
    
        # Define the filename and the extension
        filename = file.split('/')[-1].split('.')[0]
        extension = file.split('/')[-1].split('.')[1]
    
        # Tells us the trap ID and the source (filename)
        print(count + 1, filename + '.' + extension, sep=': ')
        
        # Import the .dat file
        raw_lineage = pd.read_csv(file, delimiter=' ')  # names=['time', 'division_flag', 'length', 'width', 'area', 'yfp_intensity', 'CMx', 'CMy']
    
        # Sometimes the time indices are called time and others called index
        if 'time' in raw_lineage.columns:
            raw_lineage = raw_lineage.rename(columns={'division': 'division_flag'})[['time', 'division_flag', 'length']]
        elif 'index' in raw_lineage.columns:
            raw_lineage = raw_lineage.rename(columns={'index': 'time', 'division': 'division_flag'})[['time', 'division_flag', 'length']]
    
        # So that we can correctly put in the time
        if len(np.unique(np.diff(raw_lineage['time']))) != 1 or np.unique(np.diff(raw_lineage['time']))[0] != 1:
            # print(np.unique(np.diff(raw_lineage['time'])))
            # raise IOError('time given in Wang data is not monotonous and increasing.')
            print('the time given in the data file is not increasing always by 1 so we do not know how to measure time for this lineage we will not use it.')
            continue
    
        raw_lineage['time'] = (raw_lineage['time']) / 60  # Because we map the index to the correct time-step-size which is 1 minute
    
        step_size = 1 / 60  # one-minute measurements!
    
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time']  # .round(2)
        raw_lineage['filename'] = filename
        raw_lineage = raw_lineage.reset_index(drop=True)
    
        # Add the trap ID
        raw_lineage['lineage_ID'] = count + 1 - offset
    
        if not all(x < y for x, y in zip(raw_lineage['time'].values[:-1], raw_lineage['time'].values[1:])):
            print(filename, ': Time is going backwards. We cannot use this data.')
        
            # reset the lineage_ID
            offset += 1
            continue
            
        zero_indices = raw_lineage[raw_lineage['length'] <= 0].index
        print(zero_indices)
        
        
        # irregular_step_sizes =
    
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        # If not, do not use the trace (Too much of a headache for the phenotypic variable linear regression).
        step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        if not np.all(step_sizes == step_sizes[0]):
            print(filename, ' has steps that are not the same size ', np.unique(step_sizes), ' and we are not using this data then')  # , are we are not gonna use these...')
            continue
            # exit()
    
        # Make sure there are no NaNs. If so, stop the program, something is wrong.
        if raw_lineage.isna().values.any():
            print(raw_lineage.isna().values.any())
            print(raw_lineage.isna().sum())
            print(raw_lineage)
            print(raw_lineage[raw_lineage.isnull()])
            raise IOError('there are NaNs!')


""" Create the csv files for physical, trace-centered, and trap-centered units for MM and SM data """


def main():
    import argparse
    
    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')
    
    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    # Do all the Mother and Sister Machine data
    for data_origin in wang_datasets:  # input_args.dataset_names:
        print(data_origin)
        
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
        
        # # Make sure the folders where we place the data are created already
        # create_folder(args['raw_data'])
        # create_folder(args['processed_data'])
        
        # compare_cycle_variables_to_raw_data(args)
        
        fix(args)
        
        print('*' * 200)


if __name__ == '__main__':
    main()
