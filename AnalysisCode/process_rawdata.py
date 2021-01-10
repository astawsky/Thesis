#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from AnalysisCode.global_variables import phenotypic_variables, create_folder, dataset_names, sm_datasets, mm_datasets, wang_datasets, tanouchi_datasets

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
        domain = (raw_lineage['time'].iloc[start: end].copy().values - raw_lineage['time'].iloc[start]).reshape(-1, 1)
        
        # if ppc != len(domain):
        #     print(ppc, len(domain))
        #     raise IOError('points per cycle and length of domain are not the same! {} != {}'.format(ppc, len(domain)))
        
        # domain = np.linspace(raw_lineage['time'].iloc[start], raw_lineage['time'].iloc[end], num=end - start + 1).reshape(-1, 1)
        range = np.log(raw_lineage['length'].iloc[start:end].values).reshape(-1, 1)  # the end+1 is due to indexing
        
        # Make sure they are the same size for the regression!
        assert len(domain) == len(range)
        
        if (raw_lineage['length'].iloc[start:end].values <= 0).any():
            print(raw_lineage['length'].iloc[start:end])
            plt.plot(np.arange(start, end), raw_lineage['length'].iloc[start:end], ls='--')
            plt.plot(raw_lineage['length'])
            plt.show()
            plt.close()
        
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
        
        # Check the regression
        # print(domain)
        # print(raw_lineage['time'].iloc[start:end+1].values)
        # plt.plot(domain, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.plot(domain, np.exp(range))
        
        # plt.plot(raw_lineage['time'].values, raw_lineage['length'].values, color='orange')
        # plt.scatter(raw_lineage['time'].iloc[start_indices].values, raw_lineage['length'].iloc[start_indices].values, color='blue')
        # plt.plot(raw_lineage['time'].iloc[start:end].values, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.scatter(raw_lineage['time'].iloc[start], cycle['length_birth'])
        # plt.yscale('log')
        # plt.show()
        # plt.close()
    
    # The experimental size variables
    cycle_variables_lineage['div_then_fold'] = np.append(np.nan, cycle_variables_lineage['division_ratio'].values[:-1] * np.exp(cycle_variables_lineage['fold_growth'].values[1:]))
    cycle_variables_lineage['div_and_fold'] = cycle_variables_lineage['division_ratio'] * cycle_variables_lineage['fold_growth']
    cycle_variables_lineage['fold_then_div'] = np.append(cycle_variables_lineage['division_ratio'].values[1:] * np.exp(cycle_variables_lineage['fold_growth'].values[:-1]), np.nan)
    
    # Without throwing away outliers
    without_nans = cycle_variables_lineage.copy().reset_index(drop=True)
    
    # Throwing away outliers
    cycle_variables_lineage[phenotypic_variables] = cycle_variables_lineage[phenotypic_variables].where(
        cycle_variables_lineage[phenotypic_variables] < ((3 * cycle_variables_lineage[phenotypic_variables].std()) + cycle_variables_lineage[phenotypic_variables].mean()),
        other=np.nan
    )
    
    # make them integers
    cycle_variables_lineage['lineage_ID'] = int(lin_id)
    cycle_variables_lineage['generation'] = np.arange(len(cycle_variables_lineage), dtype=int)
    cycle_variables_lineage = cycle_variables_lineage.sort_values('generation')
    
    return [cycle_variables_lineage.reset_index(drop=True), without_nans]


""" Create the csv files for physical, trace-centered, and trap-centered units for MM data of a certain type """


def main_mm(args):
    print('type of MM data we are processing:', args['data_origin'])
    
    # We have to specify the data type, which is different for each experiment
    if args['data_origin'] == 'MG1655_inLB_LongTraces':
        infiles = glob.glob(args['raw_data'] + '/*.txt')
    elif args['data_origin'] == 'Maryam_LongTraces':
        infiles = glob.glob(args['raw_data'] + '/*.csv')
        infiles = infiles + glob.glob(args['raw_data'] + '/*.xls')
    elif args['data_origin'] == 'LB_pooled':
        infiles = glob.glob(args['raw_data'] + '/*')
    elif args['data_origin'] in wang_datasets:
        infiles = glob.glob(args['raw_data'] + '/*.dat')
    
    # the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    with_outliers_cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    # The dataframe for our raw data lineages
    raw_data = pd.DataFrame(columns=['time', 'length', 'lineage_ID', 'filename'])
    
    extra_column = [
        'pos0-1',
        'pos0-1-daughter',
        'pos1',
        'pos4',
        'pos5-lower cell',
        'pos5-upper cell',
        'pos6-1-1',
        'pos6-1',
        'pos6-2',
        'pos7-1-1',
        'pos7-1-2',
        'pos7-2-1',
        'pos7-2-2',
        'pos7-3',
        'pos8',
        'pos10-1',
        'pos16-1',
        'pos16-2',
        'pos16-3',
        'pos17-1',
        'pos17-2',
        'pos17-3',
        'pos18-2',
        'pos18-3',
        'pos19-2',
        'pos19-3',
        'pos20'
    ]
    
    # In case we can't use some files we want the lineage IDs to be in integer order
    offset = 0
    
    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(infiles):
        
        filename = file.split('/')[-1].split('.')[0]
        extension = file.split('/')[-1].split('.')[1]
        
        # # Tells us the trap ID and the source (filename)
        # print(count + 1, filename + '.' + extension, sep=': ')
        
        # # creates a dataframe from the .txt or .csv file
        # if args['data_origin'] == 'MG1655_inLB_LongTraces':
        #     # The only file that has an extra column
        #     if filename.split('/')[-1] == 'pos4-4.txt':
        #         print('In this particular case the lineage divides after the first time-point and it has an extra column.')
        #         # This is because in this particular case the lineage divides after the first time-point and it has an extra column
        #         raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        #         raw_lineage = raw_lineage.iloc[1:]
        #     # All the rest don't have another column
        #     else:
        #         raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        # elif args['data_origin'] == 'Maryam_LongTraces':
        #     # This is because some the data is in .xls format while others are in .csv
        #     try:
        #         raw_lineage = pd.read_csv(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        #     except:
        #         raw_lineage = pd.read_excel(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        
        # elif args['data_origin'] == 'LAC_M9':
        #     # Simple :)
        #     raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        
        if args['data_origin'] in tanouchi_datasets:
            raw_lineage = pd.read_csv(file, delimiter=',', names=['time', 'division_flag', 'length', 'fluor', 'avg_fluor'])
            raw_lineage['time'] = (raw_lineage['time'] - 1) / 60  # Because we map the index to the correct time-step-size which is 1 minute
            step_size = 1 / 60  # one-minute measurements!
        elif args['data_origin'] == 'lambda_LB':
            # There are quite a lot of files with an extra column at the beginning
            if filename in extra_column:
                raw_lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            elif filename == 'pos15':
                print('This one is special.')
                raw_lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein', '__', '___', '___1'])[
                    ['time', 'length']]
            else:
                raw_lineage = pd.read_csv(file, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args['data_origin'] == '8-31-16 Continue':
            raw_lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args['data_origin'] in wang_datasets:
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

            raw_lineage['length'] = (raw_lineage['length']) * 0.0645  # Convert it from pixel length to micrometers
            
            step_size = 1 / 60  # one-minute measurements!
        else:
            raise IOError('This code is not meant to run the data inputted. Please label the data and put it in as an if-statement.')
        
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time']  # .round(2)
        raw_lineage['filename'] = filename
        raw_lineage = raw_lineage.reset_index(drop=True)
        # print(len(raw_lineage))
        raw_lineage = raw_lineage[raw_lineage['length'] > 0]  # Goes without saying but there was an instance of this in the Wang Data
        # print(len(raw_lineage))
        # exit()
        
        # Add the trap ID
        raw_lineage['lineage_ID'] = count + 1 - offset
        
        if not all(x < y for x, y in zip(raw_lineage['time'].values[:-1], raw_lineage['time'].values[1:])):
            print(filename, ': Time is going backwards. We cannot use this data.')
            
            # reset the lineage_ID
            offset += 1
            continue

        # if len(raw_lineage[raw_lineage['length'] <= 0]) > 0:
        #     print('This lineage has a length lower than or equal to 0, so we will not use it.')
        #     continue
        
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        # If not, do not use the trace (Too much of a headache for the phenotypic variable linear regression).
        step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        # if not np.all(step_sizes == step_sizes[0]):
        #     print(filename, ' has steps that are not the same size ', np.unique(step_sizes), ' and we are not using this data then')  # , are we are not gonna use these...')
        #     continue
        #     # exit()
        
        # This is due to data singularities
        if args['data_origin'] == '8-31-16 Continue' and .05 == step_sizes[0]:
            raw_lineage['time'] = raw_lineage['time'] * 2
            step_size = .1
        if args['data_origin'] == 'lambda_LB' and .05 == step_sizes[0]:
            pass
            step_size = .05
        elif args['data_origin'] == 'lambda_LB' and .1 == step_sizes[0]:
            print('step size == .1', filename, np.unique(step_sizes))
            raw_lineage['time'] = raw_lineage['time'] / 2
            step_size = .05
        elif args['data_origin'] == 'lambda_LB':
            print('step size == .06 repeating', filename, np.unique(step_sizes))
            raw_lineage['time'] = (raw_lineage['time'] * 15) / 20
            step_size = .05
        
        # Make sure there are no NaNs. If so, stop the program, something is wrong.
        if raw_lineage.isna().values.any():
            print(raw_lineage.isna().values.any())
            print(raw_lineage.isna().sum())
            print(raw_lineage)
            print(raw_lineage[raw_lineage.isnull()])
            raise IOError('there are NaNs!')
        
        # add it to the total data
        raw_data = raw_data.append(raw_lineage, ignore_index=True)

        # # Check division times
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'])
        # # plt.yscale('log')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # These datasets have a division flag
        if args['data_origin'] in tanouchi_datasets:  #  + wang_datasets
            start_indices = np.array(raw_lineage[raw_lineage['division_flag'] == 1].index)
            end_indices = np.append(start_indices[1:] - 1, len(raw_lineage) - 1)
            # exit()
        else:
            # Figure out the indices for the division events
            start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # the inter-division times
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        # Each cycle must be longer than 10 minutes as a physical constraint to help us with the division events
        at_least_number = np.where(cycle_durations > .1)[0]
        start_indices, end_indices = start_indices[at_least_number], end_indices[at_least_number]
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        if len(cycle_durations) != len(at_least_number):
            print('{} cycle(s) were faster than 10 minutes, so we will not use them'.format(len(cycle_durations) - len(at_least_number)))
        
        # # Apply the four data points condition
        # start_indices, end_indices, cycle_durations = start_indices[at_least_number], end_indices[at_least_number], cycle_durations[at_least_number]
        
        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / step_size), dtype=int)
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage, with_outliers = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, int(count + 1 - offset), fit_the_lengths=True)

        # # Check division times
        # interpolation = np.array([lb * np.exp(np.linspace(0, gt, int(np.round(gt * 60, 0)) + 1) * gr) for gt, gr, lb in zip(with_outliers['generationtime'], with_outliers['growth_rate'], with_outliers['length_birth'])])
        # interpolation = np.concatenate(interpolation, axis=0)
        #
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'], color='blue', label='raw_lineage')
        # plt.plot(np.arange(start_indices[0], len(interpolation)+start_indices[0]), interpolation, color='orange', label='cycle variable fit', ls='--')
        # plt.scatter(start_indices, raw_lineage['length'].iloc[start_indices], color='green', label='start indices')
        # plt.scatter(end_indices, raw_lineage['length'].iloc[end_indices], color='red', label='end indices')
        # plt.title(filename)
        # plt.xlabel('index')
        # plt.ylabel(r'length ($\mu$m)')
        # plt.yscale('log')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # append the cycle variables to the
        cycle_variables = cycle_variables.append(cycle_variables_lineage, ignore_index=True)
        with_outliers_cycle_variables = with_outliers_cycle_variables.append(with_outliers, ignore_index=True)
    
    print('processed data:\n', cycle_variables)
    print('cleaned raw data:\n', raw_data)
    
    for variable in phenotypic_variables:
        if variable == 'division_ratio':
            print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - len(cycle_variables['lineage_ID'].unique())))
        else:
            print(variable, cycle_variables[variable].count() / len(cycle_variables[variable]))
    
    # reset the index for good practice
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(os.path.dirname(os.path.dirname(args['raw_data'])) + '/raw_data_all_in_one.csv', index=False)
    # raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(args['processed_data'] + '/raw_data.csv', index=False)
    # minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['lineage_ID']).to_csv(args['processed_data'] + '/raw_data_tc.csv', index=False)
    
    without_outliers = args['processed_data'] + 'z_score_under_3'
    create_folder(without_outliers)
    
    cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(without_outliers + '/physical_units_without_outliers.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        without_outliers + '/trace_centered_without_outliers.csv',
        index=False)
    
    with_outliers_cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args['processed_data'] + '/physical_units.csv', index=False)
    minusing(with_outliers_cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        args['processed_data'] + '/trace_centered.csv',
        index=False)


""" Recreate the raw data from smoothed exponential regression and see how well they compare """


def compare_cycle_variables_to_raw_data(args):
    import seaborn as sns
    
    print(args['data_origin'])
    
    # import the labeled measured bacteria in physical units and the raw data in physical units
    physical_units = pd.read_csv('{}/physical_units.csv'.format(args['processed_data']))
    raw_data = pd.read_csv('{}/raw_data.csv'.format(args['processed_data']))
    
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
            # plt.title(args['data_origin'] + ': ' + str(lin_id))
            # plt.show()
            # plt.close()
        
        # recreated = [x_0 * np.exp(alpha * time) for time in np.arange(0, tau + .05, .05) for tau, alpha, x_0 in zip(lineage['generationtime'], lineage['growth_rate'], lineage['length_birth'])]
        
        x = raw_trace['time'][start_indices[0]:end_indices[-1] + 1]
        # x = np.linspace(raw_trace['time'].values[start], raw_trace['time'].values[end], num=end - start + 1)
        
        # print(len(x), len(raw_trace['length'][start_indices[0]:end_indices[-1] + 1]))
        # print(len(recreated))
        # print('----')
        # print(x, raw_trace['length'][start_indices[0]:end_indices[-1] + 1], sep='\n')
        # exit()
        
        plt.plot(x, raw_trace['length'][start_indices[0]:end_indices[-1] + 1], label='raw data')
        # plt.plot(x, recreated, label='recreated from regression')
        plt.legend()
        plt.tight_layout()
        plt.title(args['data_origin'] + ': ' + str(lin_id))
        plt.show()
        plt.close()


""" Create the csv files for physical, trace-centered, and trap-centered units for SM data """


def main_sm(args):
    # Where we will put the raw data
    raw_data = pd.DataFrame(columns=['time', 'length', 'dataset', 'trap_ID', 'trace', 'lineage_ID'])
    
    # The files that have the RawData
    files = glob.glob(args['raw_data'] + '/*.xls')
    
    # load first sheet of each Excel-File, fill rawdata dataframe
    for count, file in enumerate(files):
        # print(count, file.split('/')[-1], sep=': ')
        
        # print(45 % 5)
        # print(46 % 5)
        # exit()
        
        # creates a dataframe from the excel file
        tmpdata = pd.read_excel(file)
        
        # Make sure there are no NaNs in the data
        assert ~tmpdata.isna().values.any()
        
        # Determine the relationship between the two lineage in a file depending on the name of the file
        if ('sis' in file.split('/')[-1]) or ('SIS' in file.split('/')[-1]):
            dataset = 'SL'
        else:
            dataset = 'NL'
        
        # Separate and categorize the traces in the trap
        if dataset == 'SL':
            a_trace = tmpdata[['timeA', 'lengthA']].rename(columns={'timeA': 'time', 'lengthA': 'length'})
            b_trace = tmpdata[['timeB', 'lengthB']].rename(columns={'timeB': 'time', 'lengthB': 'length'})
        else:
            a_trace = tmpdata[['timeA', 'L1']].rename(columns={'timeA': 'time', 'L1': 'length'})
            b_trace = tmpdata[['timeB', 'L2']].rename(columns={'timeB': 'time', 'L2': 'length'})
        
        # Are they SL or NL?
        a_trace['dataset'] = dataset
        b_trace['dataset'] = dataset
        
        # What trap are they in from the pooled ensemble?
        a_trace['trap_ID'] = (count + 1)
        b_trace['trap_ID'] = (count + 1)
        
        # Arbitrarily name the traces
        a_trace['trace'] = 'A'
        b_trace['trace'] = 'B'
        
        # Give each lineage a unique ID
        a_trace['lineage_ID'] = (count + 1) * 2 - 1
        b_trace['lineage_ID'] = (count + 1) * 2
        
        # Set the floats to be accurate to the 2nd decimal point because of the timesteps in hours
        a_trace['time'] = a_trace['time'].round(2)
        b_trace['time'] = b_trace['time'].round(2)
        
        # Check that they didn't round to 0.06 or something
        if any([int(np.round(l * 100, 0) % 5) for l in a_trace['time'].values]) or any([int(np.round(l * 100, 0) % 5) for l in a_trace['time'].values]):
            print('Rounded Wrong to not .05 multiples')
            exit()
        
        # Check if time is going forward for the "A" trace
        time_monotony_a = all(x < y for x, y in zip(a_trace['time'].values[:-1], a_trace['time'].values[1:]))
        # Check if time is going forward for the "B" trace
        time_monotony_b = all(x < y for x, y in zip(b_trace['time'].values[:-1], b_trace['time'].values[1:]))
        
        if (not time_monotony_a) or (not time_monotony_b):
            print(file, ': Time is going backwards. We cannot use this data.')
            print("False is bad! --", "A:", time_monotony_a, "B:", time_monotony_b)
            
            # # To see this in the lineage itself
            # plt.plot(raw_lineage['time'], raw_lineage['length'])
            # plt.show()
            # plt.close()
            
            continue
        
        # the data contains all dataframes from the excel files in the directory _infiles
        raw_data = raw_data.append(a_trace, ignore_index=True)
        raw_data = raw_data.append(b_trace, ignore_index=True)
    
    # There must be some data
    assert len(raw_data) > 0
    # There can't be any NaNs
    assert ~raw_data.isna().values.any()
    
    # # Save the raw data to .csv format
    # raw_data.sort_values(['dataset', 'trap_ID', 'trace']).reset_index(drop=True).to_csv(args['processed_data'] + '/raw_data.csv', index=False)
    # minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['dataset', 'trap_ID', 'trace']).reset_index(drop=True).to_csv(args['processed_data'] + '/raw_data_tc.csv', index=False)
    
    # This is the order of the cycle variables in the processed dataframe
    order = phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'lineage_ID', 'generation']
    
    # The dataframe for our variables
    cycle_variables = pd.DataFrame(columns=order)
    with_outliers_cycle_variables = pd.DataFrame(columns=order)
    
    for lineage_id in raw_data.lineage_ID.unique():
        # print('Lineage ID:', lineage_id)
        
        # Get the lineage
        raw_lineage = raw_data[raw_data['lineage_ID'] == lineage_id]
        
        # Figure out the indices for the division events
        start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # the inter-division times
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        # Each cycle must consist of at least four data points
        at_least_number = np.where(cycle_durations > (2 * 0.05))[0]
        
        # Apply the four data points condition
        start_indices, end_indices, cycle_durations = start_indices[at_least_number], end_indices[at_least_number], cycle_durations[at_least_number]
        
        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / .05) + np.ones_like(cycle_durations), dtype=int)
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage, with_outliers = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, int(lineage_id), fit_the_lengths=True)

        # # Check division times
        # interpolation = np.array([lb * np.exp(np.linspace(0, gt, int(np.round(gt * 20, 0)) + 1) * gr) for gt, gr, lb in zip(with_outliers['generationtime'], with_outliers['growth_rate'], with_outliers['length_birth'])])
        # interpolation = np.concatenate(interpolation, axis=0)
        #
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'].values, color='blue', label='raw_lineage')
        # plt.plot(np.arange(start_indices[0], len(interpolation)+start_indices[0]), interpolation, color='orange', label='cycle variable fit', ls='--')
        # plt.scatter(start_indices, raw_lineage['length'].iloc[start_indices], color='green', label='start indices')
        # plt.scatter(end_indices, raw_lineage['length'].iloc[end_indices], color='red', label='end indices')
        # plt.title(lineage_id)
        # plt.xlabel('index')
        # plt.ylabel(r'length ($\mu$m)')
        # plt.yscale('log')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # Add the SM categorical variables
        cycle_variables_lineage['trap_ID'] = raw_lineage['trap_ID'].unique()[0]
        cycle_variables_lineage['trace'] = raw_lineage['trace'].unique()[0]
        cycle_variables_lineage['dataset'] = raw_lineage['dataset'].unique()[0]
        
        with_outliers['trap_ID'] = raw_lineage['trap_ID'].unique()[0]
        with_outliers['trace'] = raw_lineage['trace'].unique()[0]
        with_outliers['dataset'] = raw_lineage['dataset'].unique()[0]
        
        # Append the cycle variables to the processed dataframe
        cycle_variables = cycle_variables.append(cycle_variables_lineage[order], ignore_index=True)
        with_outliers_cycle_variables = with_outliers_cycle_variables.append(with_outliers[order], ignore_index=True)
    
    print('processed data:\n', cycle_variables)
    
    print('cleaned raw data:\n', raw_data)
    
    # Check how much NaNs were introduced because of the zscore < 3 condition on one of the dataframes (no outliers)
    for variable in phenotypic_variables:
        if variable in ['division_ratio', 'div_and_fold', 'fold_then_div']:
            # This is because, by definition, the first two variables have 1 NaN value at the first generation, while the third variable has 1 NaN at the end of each lineage
            print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - (1 * len(cycle_variables['lineage_ID'].unique()))))
        elif variable in ['div_then_fold']:
            # This is because, by definition, div then fold variable has two NaNs at the first two generations of each lineage
            print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - (2 * len(cycle_variables['lineage_ID'].unique()))))
        else:
            print(variable, cycle_variables[variable].count() / len(cycle_variables[variable]))
    
    # reset the index for good practice
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(os.path.dirname(os.path.dirname(args['raw_data'])) + '/raw_data_all_in_one.csv', index=False)
    
    without_outliers = args['processed_data'] + 'z_score_under_3'
    create_folder(without_outliers)
    
    cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(without_outliers + '/physical_units_without_outliers.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        without_outliers + '/trace_centered_without_outliers.csv',
        index=False)
    
    with_outliers_cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args['processed_data'] + '/physical_units.csv', index=False)
    minusing(with_outliers_cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        args['processed_data'] + '/trace_centered.csv',
        index=False)


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
        
        if data_origin in sm_datasets:  # Get SM data
            # Get SM data
            main_sm(args)
        else:
            # Get MM data
            main_mm(args)
        
        print('*' * 200)
        exit()
        

if __name__ == '__main__':
    main()
