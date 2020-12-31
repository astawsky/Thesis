#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, dataset_names

""" This function takes care of which statistic we want to subtract the trajectories to get what we asked for """


def minusing(info, parameters):
    tc = info.copy()
    for lineage_id in info['lineage_ID'].unique():
        trace = info[(info['lineage_ID'] == lineage_id)].copy()
        tc.loc[trace.index, parameters] = trace[parameters] - trace[parameters].mean()
    
    return tc


""" Recognizes where division events occur """


def get_division_indices(raw_trace):
    # From the raw data, we see where the difference between two values of length
    # falls drastically, suggesting a division has occurred.
    diffs = np.diff(raw_trace)
    
    # How much of a difference there has to be between two points to consider division having taken place.
    threshold_of_difference_for_division = -1.2
    
    # The array of indices where division took place.
    index_div = np.where(diffs < threshold_of_difference_for_division)[0].flatten()
    
    # If the two consecutive indices in the array of indices where division took place are
    # less than two time steps away, then we discard them
    for ind1, ind2 in zip(index_div[:-1], index_div[1:]):
        if ind2 - ind1 <= 2:
            index_div = np.delete(index_div, np.where(index_div == ind2))
    
    # An array of indices where a cycle starts and ends
    start_indices = [x + 1 for x in index_div]
    end_indices = [x for x in index_div]
    
    # Make sure they are the same size
    assert len(start_indices) == len(end_indices)
    
    # if raw_trace[0] in start_indices:
    #     del start_indices[np.where(start_indices != raw_trace[0])[0]]
    
    # Count the first cycle by adding 0
    if 0 not in start_indices:
        start_indices.append(0)  # to count the first cycle
        start_indices.sort()  # because '0' above was added at the end and the rest are in order already
    del start_indices[-1]  # we do not count the last cycle because it most likely did not finish by the time recording ended
    end_indices.sort()  # same as with start_indices, just in case
    
    # If the last index is in the array of indices where a cycle starts, then it obviously a mistake and must be removed
    if raw_trace[-1] in start_indices:
        start_indices.remove(raw_trace[-1])
    if 1 in start_indices and 0 in start_indices:
        start_indices.remove(0)
    
    # Similarly, if the fist index '0' is in the array of indices where a cycle ends, then it obviously is a mistake and must be removed
    if 0 in end_indices:
        end_indices.remove(0)
    
    # Sanity check: If the starting index for a cycle comes after the ending index, then it is obviously a mistake
    # Not a foolproof method so, if the length of the bacteria drops after the start time of the cycle, we know it is not the real starting
    # time since a bacteria does not shrink after dividing
    for start, end in zip(start_indices, end_indices):
        if start >= end:
            IOError('start', start, 'end', end, "didn't work")
    
    # Make them numpy arrays now so we can subtract them
    start_indices = np.array(start_indices)
    end_indices = np.array(end_indices)
    
    if len(start_indices) != len(end_indices):
        print(start_indices)
        print(end_indices)
        print(len(start_indices))
        print(len(end_indices))
        exit()
    
    # Make sure they are the same size
    assert len(start_indices) == len(end_indices)
    
    return [start_indices, end_indices]


"""We use this for linear regression in the cycle parameter process """


def linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, lin_id, fit_the_lengths):
    # Make sure everything is aligned
    assert len(start_indices) == len(cycle_durations) == len(data_points_per_cycle) == len(end_indices)
    
    # the dataframe for our variables
    cycle_variables_lineage = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    for start, end, inter_div_time, gen in zip(start_indices, end_indices, cycle_durations, np.arange(len(cycle_durations), dtype=int)):
        
        # for our regression
        domain = (raw_lineage['time'].iloc[start: end + 1].copy().values - raw_lineage['time'].iloc[start]).reshape(-1, 1)
        
        # domain = np.linspace(raw_lineage['time'].iloc[start], raw_lineage['time'].iloc[end], num=end - start + 1).reshape(-1, 1)
        range = np.log(raw_lineage['length'].iloc[start:end + 1].values).reshape(-1, 1)  # the end+1 is due to indexing
        
        # Make sure they are the same size for the regression!
        assert len(domain) == len(range)
        
        # our regression
        reg = LinearRegression().fit(domain, range)
        
        # the phenotypic variables of a cycle
        cycle = pd.Series()
        
        # Define the cycle phenotypic variables
        cycle['generationtime'] = inter_div_time
        cycle['growth_rate'] = reg.coef_[0][0]
        cycle['fold_growth'] = cycle['growth_rate'] * cycle['generationtime']
        
        # Categorical labels to identify lineage and a cycle in it
        cycle['lineage_ID'] = int(lin_id)
        cycle['generation'] = gen
        
        # Do the length at birth and length at division come straight from the data or from the regression?
        if fit_the_lengths:
            # phenotypic variables
            cycle['length_birth'] = np.exp(reg.predict(domain[0].reshape(-1, 1))[0][0])
            cycle['length_final'] = np.exp(reg.predict(domain[-1].reshape(-1, 1))[0][0])
        else:
            # phenotypic variables
            cycle['length_birth'] = np.exp(range[0][0])
            cycle['length_final'] = np.exp(range[-1][0])
        
        # We define the division ratio as how much length a cell received from its mother,
        # since the first generation has no recorded mother the value will be a NaN.
        if len(cycle_variables_lineage) == 0:
            cycle['division_ratio'] = np.nan
        else:
            cycle['division_ratio'] = cycle['length_birth'] / cycle_variables_lineage['length_final'].iloc[-1]
            # try:
            #     cycle['division_ratio'] = cycle['length_birth'] / cycle_variables_lineage['length_final'].iloc[-1]
            # except:
            #     print('+'*200, start, cycle['length_birth'], cycle_variables_lineage, sep='\n')
            #     exit()
        
        # After defining the lengths, get the added length
        cycle['added_length'] = cycle['length_final'] - cycle['length_birth']
        
        # Add them to the cycle variables of the lineage
        cycle_variables_lineage = cycle_variables_lineage.append(cycle, ignore_index=True)
        
        # print(start)
        # print(cycle)
        # print(cycle_variables_lineage)
        
        # # Check the regression
        # plt.plot(domain, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.plot(domain, np.exp(range))
        # plt.show()
        # plt.close()
    
    # Without throwing away outliers
    without_nans = cycle_variables_lineage.copy().reset_index(drop=True)
    
    # Throwing away outliers
    cycle_variables_lineage[phenotypic_variables] = cycle_variables_lineage[phenotypic_variables].where(
        cycle_variables_lineage[phenotypic_variables] < ((3 * cycle_variables_lineage[phenotypic_variables].std()) + cycle_variables_lineage[phenotypic_variables].mean()),
        other=np.nan
    )
    
    # Categorical variables, make them integers
    cycle_variables_lineage['lineage_ID'] = int(lin_id)
    cycle_variables_lineage['generation'] = np.arange(len(cycle_variables_lineage), dtype=int)
    cycle_variables_lineage = cycle_variables_lineage.sort_values('generation')
    
    return [cycle_variables_lineage.reset_index(drop=True), without_nans]


""" Create the csv files for physical, trace-centered, and trap-centered units for MM data of a certain type """


def main_mm(args):
    print('type of MM data we are processing:', args.data_origin)
    
    # We have to specify the data type, which is different for each experiment
    if args.data_origin == 'MG1655_inLB_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.txt')
    elif args.data_origin == 'Maryam_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.csv')
        infiles = infiles + glob.glob(args.raw_data + '/*.xls')
    elif args.data_origin == 'LB_pooled':
        infiles = glob.glob(args.raw_data + '/*')
    else:
        infiles = glob.glob(args.raw_data + '/*.txt')
    
    # the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    with_outliers_cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    # The dataframe for our raw data lineages
    raw_data = pd.DataFrame(columns=['time', 'length', 'lineage_ID', 'filename'])
    
    extra_column = [
        'pos0-1.txt',
        'pos0-1-daughter.txt',
        'pos1.txt',
        'pos4.txt',
        'pos5-lower cell.txt',
        'pos5-upper cell.txt',
        'pos6-1-1.txt',
        'pos6-1.txt',
        'pos6-2.txt',
        'pos7-1-1.txt',
        'pos7-1-2.txt',
        'pos7-2-1.txt',
        'pos7-2-2.txt',
        'pos7-3.txt',
        'pos8.txt',
        'pos10-1.txt',
        'pos16-1.txt',
        'pos16-2.txt',
        'pos16-3.txt',
        'pos17-1.txt',
        'pos17-2.txt',
        'pos17-3.txt',
        'pos18-2.txt',
        'pos18-3.txt',
        'pos19-2.txt',
        'pos19-3.txt',
        'pos20.txt'
    ]
    
    # When we take out the outliers we keep trace of how many cycles we lost wrt. the lineage and the pooled ensemble
    whats_left = {'variable': [], 'lineage': [], 'number_taken_out': []}
    
    # In case we can't use some files we want the lineage IDs to be in integer order
    offset = 0
    
    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(infiles):
        
        # So we know how to read it
        filename = file.split('/')[-1].split('.')[0]
        extension = file.split('/')[-1].split('.')[1]
        
        # Tells us the trap ID and the source (filename)
        print(count, filename + '.' + extension, sep=': ')
        
        # creates a dataframe from the .txt or .csv file
        # if args.data_origin == 'MG1655_inLB_LongTraces':
        #     # The only file that has an extra column
        #     if filename.split('/')[-1] == 'pos4-4.txt':
        #         print('In this particular case the lineage divides after the first time-point and it has an extra column.')
        #         # This is because in this particular case the lineage divides after the first time-point and it has an extra column
        #         raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        #         raw_lineage = raw_lineage.iloc[1:]
        #     # All the rest don't have another column
        #     else:
        #         raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        # elif args.data_origin == 'Maryam_LongTraces':
        #     # This is because some the data is in .xls format while others are in .csv
        #     try:
        #         raw_lineage = pd.read_csv(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        #     except:
        #         raw_lineage = pd.read_excel(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        # if args.data_origin == 'lambda_LB':
        #     # There are quite a lot of files with an extra column at the beginning
        #     if filename in extra_column:
        #         raw_lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        #     elif filename.split('/')[-1] == 'pos15.txt':
        #         print('This one is special.')
        #         raw_lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein', '__', '___', '___1'])[
        #             ['time', 'length']]
        #     else:
        #         raw_lineage = pd.read_csv(file, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        # elif args.data_origin == 'LAC_M9':
        #     # Simple :)
        #     raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        if args.data_origin in ['MC4100_25C', 'MC4100_27C', 'MC4100_37C']:
            raw_lineage = pd.read_csv(file, delimiter=',', names=['time', 'division_flag', 'length', 'fluor', 'avg_fluor'])
            raw_lineage['time'] = (raw_lineage['time']-1) / 60
        else:
            raise IOError('This code is not meant to run the data inputted. Please label the data and put it in as an if-statement.')
        
        print(raw_lineage)
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time']
        raw_lineage['filename'] = filename
        raw_lineage = raw_lineage.reset_index(drop=True)
        
        # Add the trap ID
        raw_lineage['lineage_ID'] = count - offset
        
        if not all(x < y for x, y in zip(raw_lineage['time'].values[:-1], raw_lineage['time'].values[1:])):
            print(filename, ': Time is going backwards. We cannot use this data.')
            
            # reset the lineage_ID
            offset += 1
            
            # for x, y in zip(raw_lineage['time'].values[:-1], raw_lineage['time'].values[1:]):
            #     if x >= y:
            #         print(x, y)
            #
            # print('~' * 200)
            #
            # # To see this in the lineage itself
            # plt.plot(raw_lineage['time'], raw_lineage['length'])
            # plt.show()
            # plt.close()
            
            continue
        
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        # If not, do not use the trace (Too much of a headache for the phenotypic variable linear regression).
        step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        if not np.all(step_sizes == step_sizes[0]):
            print(filename, ' has steps that are not the same size')  # , are we are not gonna use these...')
            # continue
        
        # Make sure there are no NaNs. If so, stop the program, something is wrong.
        if raw_lineage.isna().values.any():
            print(raw_lineage.isna().values.any())
            print(raw_lineage.isna().sum())
            print(raw_lineage)
            print(raw_lineage[raw_lineage.isnull()])
            exit()
        
        # print('unique step size:', np.unique(step_sizes))
        # print(raw_lineage)
        
        # # Figure out the indices for the division events
        # start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # print(len(start_indices), raw_lineage['division_flag'].sum())
        # print(start_indices, list(raw_lineage[raw_lineage['division_flag'] == 1].index))
        
        # add it to the total data
        raw_data = raw_data.append(raw_lineage, ignore_index=True)
        
        # # Check division times
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'])
        # for start, end in zip(start_indices, end_indices):
        #     plt.axvline(start, color='green')
        #     plt.axvline(end, color='red')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
    print('cleaned raw data:\n', raw_data)
    
    # reset the index for good practice
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data.csv', index=False)
    minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data_tc.csv', index=False)


""" Recreate the raw data from smoothed exponential regression and see how well they compare """


def compare_cycle_variables_to_raw_data(args):
    import seaborn as sns
    
    print(args.data_origin)
    
    # import the labeled measured bacteria in physical units and the raw data in physical units
    physical_units = pd.read_csv('{}/physical_units.csv'.format(args.save_folder))
    raw_data = pd.read_csv('{}/raw_data.csv'.format(args.save_folder))
    
    # Check that we have the same amount of lineages in each dataframe
    assert all([lp == lr for lp, lr in zip(physical_units.lineage_ID.unique(), raw_data.lineage_ID.unique())])
    
    for lin_id in physical_units.lineage_ID.unique():
        # define the lineages we are working with
        lineage = physical_units[physical_units['lineage_ID'] == lin_id]
        start_indices, end_indices = get_division_indices(raw_data[raw_data['lineage_ID'] == lin_id]['length'].values)
        
        # Compare only the part of the raw trace that we computed the cycle variables for
        raw_trace = raw_data[raw_data['lineage_ID'] == lin_id].sort_values('time')
        
        # stylistic reasons
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        recreated = []
        print(raw_trace)
        for start, end, alpha, x_0 in zip(start_indices, end_indices, lineage['growth_rate'], lineage['length_birth']):
            # cycle = []
            # print(tau, alpha, x_0)
            
            # get the generationtime
            tau = raw_trace['time'].values[end] - raw_trace['time'].values[start]
            
            for time in np.linspace(0, tau, num=end - start + 1):
                # print('time:', time)
                # print('size:', x_0 * np.exp(alpha * time))
                recreated.append(x_0 * np.exp(alpha * time))
                # cycle.append(x_0 * np.exp(alpha * time))
            
            # x = np.linspace(raw_trace['time'].values[start], raw_trace['time'].values[end], num=end-start+1)
            # print(len(x))
            # print(len(raw_trace[['time', 'length']][start:end+1].values))
            # plt.plot(x, raw_trace['length'][start:end+1].values, label='raw data')
            # plt.plot(x, cycle, label='recreated from regression')
            # plt.legend()
            # plt.tight_layout()
            # plt.title(args.data_origin + ': ' + str(lin_id))
            # plt.show()
            # plt.close()
        
        # recreated = [x_0 * np.exp(alpha * time) for time in np.arange(0, tau + .05, .05) for tau, alpha, x_0 in zip(lineage['generationtime'], lineage['growth_rate'], lineage['length_birth'])]
        
        x = raw_trace['time'][start_indices[0]:end_indices[-1] + 1]
        # x = np.linspace(raw_trace['time'].values[start], raw_trace['time'].values[end], num=end - start + 1)
        
        # print(len(x), len(raw_trace['length'][start_indices[0]:end_indices[-1] + 1]))
        # print(len(recreated))
        # print('----')
        print(x, raw_trace['length'][start_indices[0]:end_indices[-1] + 1], sep='\n')
        # exit()
        
        plt.plot(x, raw_trace['length'][start_indices[0]:end_indices[-1] + 1], label='raw data')
        # plt.plot(x, recreated, label='recreated from regression')
        plt.legend()
        plt.tight_layout()
        plt.title(args.data_origin + ': ' + str(lin_id))
        plt.show()
        plt.close()


""" Create the csv files for physical, trace-centered, and trap-centered units for MM and SM data """


def main():
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in dataset_names:
        
        # Create the arguments for this function
        parser = argparse.ArgumentParser(description='Process Mother Machine and Sister Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        
        # We need this for the SM data processing
        if data_origin == 'SM':
            parser.add_argument('-SL', '--sl_infiles', metavar='', type=str,
                                help='Location of Sister Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/SL')
            parser.add_argument('-NL', '--nl_infiles', metavar='', type=str,
                                help='Location of Neighboring Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/NL')
        
        # Finalize the arguments
        args = parser.parse_args()
        
        # Make sure the folders where we place the data are created already
        create_folder(args.raw_data)
        create_folder(args.save_folder)
        
        # compare_cycle_variables_to_raw_data(args)
        
        main_mm(args)
        
        print('*' * 200)
    
    # # How long did it take to do the mother machine?
    # mm_time = time.time() - first_time
    #
    # # Now we move onto the Sister Machine data
    # data_origin = 'SM'
    #
    # parser = argparse.ArgumentParser(description='Process Sister Machine Lineage Data.')
    # parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    # parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
    #                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin)
    # parser.add_argument('-SL', '--sl_infiles', metavar='', type=str,
    #                     help='Location of Sister Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/SL')
    # parser.add_argument('-NL', '--nl_infiles', metavar='', type=str,
    #                     help='Location of Neighboring Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/NL')
    # parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
    #                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    # args = parser.parse_args()
    #
    # create_folder(args.raw_data)
    # create_folder(args.save_folder)
    #
    # main_sm(args)
    #
    # print('*' * 200)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))


if __name__ == '__main__':
    main()
