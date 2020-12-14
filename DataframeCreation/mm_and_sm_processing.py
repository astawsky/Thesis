#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
from sklearn.linear_model import LinearRegression
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, mm_data_names

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
    
    # Make sure they are the same size
    assert len(start_indices) == len(end_indices)
    
    return [start_indices, end_indices]


"""We use this for linear regression in the cycle parameter process """


def linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, trap_id, fit_the_lengths):
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
        cycle['lineage_ID'] = trap_id
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
    
    return cycle_variables_lineage


""" Create the csv files for physical, trace-centered, and trap-centered units for MM data of a certain type """


def main_mm(args):
    print('type of MM data we are processing:', args.data_origin)
    
    # We have to specify the data type, which is different for each experiment
    if args.data_origin == 'MG1655_inLB_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.txt')
    elif args.data_origin == 'Maryam_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.csv')
        infiles = infiles + glob.glob(args.raw_data + '/*.xls')
    else:
        infiles = glob.glob(args.raw_data + '/*.txt')
    
    # the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    # The dataframe for our raw data lineages
    raw_data = pd.DataFrame(columns=['time', 'length', 'lineage_ID'])
    
    # filenames_with_nans = [
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy02ch03R.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy03ch02.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy03ch01.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy09ch05.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy08ch01.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_xy11ch01.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_xy11ch03.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_xy11ch06.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_xy11ch07.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_xy09ch02.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy02ch08L.csv',
    #     '/Users/alestawsky/PycharmProjects/Thesis/RawData/Maryam_LongTraces/d_09252019_nd041xy06ch01.csv'
    # ]
    
    extra_column = [
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos0-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos0-1-daughter.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos4.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos5-lower cell.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos5-upper cell.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos6-1-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos6-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos6-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos7-1-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos7-1-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos7-2-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos7-2-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos7-3.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos8.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos10-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos15.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos16-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos16-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos16-3.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos17-1.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos17-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos17-3.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos18-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos18-3.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos19-2.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos19-3.txt',
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos20.txt'
    ]
    
    # load first sheet of each Excel-File, fill internal data structure]
    for count, filename in enumerate(infiles):
        
        # Tells us the trap ID and the source (filename)
        print(count, filename.split('/')[-1], sep=': ')
        
        # creates a dataframe from the .txt or .csv file
        if args.data_origin == 'MG1655_inLB_LongTraces':
            # The only file that has an extra column
            if filename.split('/')[-1] == 'pos4-4.txt':
                print('In this particular case the lineage divides after the first time-point and it has an extra column.')
                # This is because in this particular case the lineage divides after the first time-point and it has an extra column
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
                raw_lineage = raw_lineage.iloc[1:]
            # All the rest don't have another column
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args.data_origin == 'Maryam_LongTraces':
            # This is because some the data is in .xls format while others are in .csv
            try:
                raw_lineage = pd.read_csv(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
            except:
                raw_lineage = pd.read_excel(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        elif args.data_origin == 'lambda_LB':
            # There are quite a lot of files with an extra column at the beginning
            if filename in extra_column:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args.data_origin == 'LAC_M9':
            # Simple :)
            raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        else:
            raise IOError('This code is not meant to run the data inputted. Please label the data and put it in as an if-statement.')
        
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time'].round(2)
        raw_lineage = raw_lineage.reset_index(drop=True)
        
        # Add the trap ID
        raw_lineage['lineage_ID'] = count
        
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        # If not, do not use the trace (Too much of a headache for the phenotypic variable linear regression).
        step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        if not np.all(step_sizes == step_sizes[0]):
            print(filename, ' has steps that are not the same size, are we are not gonna use these...')
            continue
        
        # Make sure there are no NaNs. If so, stop the program, something is wrong.
        if raw_lineage.isna().values.any():
            print(raw_lineage.isna().values.any())
            print(raw_lineage.isna().sum())
            print(raw_lineage)
            print(raw_lineage[raw_lineage.isnull()])
            exit()
        
        # add it to the total data
        raw_data = raw_data.append(raw_lineage, ignore_index=True)
        
        # Figure out the indices for the division events
        start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # # Check division times
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'])
        # for start, end in zip(start_indices, end_indices):
        #     plt.axvline(start, color='green')
        #     plt.axvline(end, color='red')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # the inter-division times
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        # Each cycle must consist of at least four data points
        at_least_number = np.where(cycle_durations > (2 * 0.05))[0]
        
        # Apply the four data points condition
        start_indices, end_indices, cycle_durations = start_indices[at_least_number], end_indices[at_least_number], cycle_durations[at_least_number]
        
        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / .05) + np.ones_like(cycle_durations), dtype=int)
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, int(count), fit_the_lengths=True)
        
        # append the cycle variables to the
        cycle_variables = cycle_variables.append(cycle_variables_lineage, ignore_index=True)
        
        # # This is to check that the point are where they need to be
        # for variable in phenotypic_variables:
        #     plt.plot(cycle_variables[cycle_variables['trap_ID'] == count][variable])
        #     plt.axhline(cycle_variables[cycle_variables['trap_ID'] == count][variable].mean())
        #     plt.title(variable)
        #     plt.show()
        #     plt.close()
    
    print('processed data:\n', cycle_variables)
    print('cleaned raw data:\n', raw_data)
    
    # reset the index for good practice
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data.csv', index=False)
    cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args.save_folder + '/physical_units.csv', index=False)
    minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data_tc.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args.save_folder + '/trace_centered.csv',
                                                                                                                                                   index=False)


""" Reading the raw data into a big pandas Dataframes """


def get_sm_rawdata(infiles, dataset, data, offset=0):
    # load first sheet of each Excel-File, fill internal data structure]
    for count, filename in enumerate(infiles):
        print(count + offset, filename.split('/')[-1], sep=': ')
        
        # creates a dataframe from the excel file
        tmpdata = pd.read_excel(filename)
        
        # Make sure there are no NaNs in the data
        assert ~tmpdata.isna().values.any()
        
        # Separate and categorize the traces in the trap
        if dataset == 'SL':
            a_trace = tmpdata[['timeA', 'lengthA']].rename(columns={'timeA': 'time', 'lengthA': 'length'})
            b_trace = tmpdata[['timeB', 'lengthB']].rename(columns={'timeB': 'time', 'lengthB': 'length'})
        else:
            a_trace = tmpdata[['timeA', 'L1']].rename(columns={'timeA': 'time', 'L1': 'length'})
            b_trace = tmpdata[['timeB', 'L2']].rename(columns={'timeB': 'time', 'L2': 'length'})
        a_trace['dataset'] = dataset
        b_trace['dataset'] = dataset
        a_trace['trap_ID'] = (count + offset + 1)
        b_trace['trap_ID'] = (count + offset + 1)
        a_trace['trace'] = 'A'
        b_trace['trace'] = 'B'
        a_trace['lineage_ID'] = (count + offset + 1) * 2 - 1
        b_trace['lineage_ID'] = (count + offset + 1) * 2
        
        # the data contains all dataframes from the excel files in the directory _infiles
        data = data.append(a_trace, ignore_index=True)
        data = data.append(b_trace, ignore_index=True)
    
    # There must be some data
    assert len(data) > 0
    # There can't be any NaNs
    assert ~data.isna().values.any()
    
    return data


""" Create the csv files for physical, trace-centered, and trap-centered units for SM data """


def main_sm(args):
    # Where we will put the raw data
    raw_data = pd.DataFrame(columns=['time', 'length', 'dataset', 'trap_ID', 'trace', 'lineage_ID'])
    
    # for dataset, infiles in zip(['SL', 'NL'], [glob.glob(args.sl_infiles + '/*.xls'), glob.glob(args.nl_infiles + '/*.xls')]):
    #     print(dataset)
    #     raw_data = get_sm_rawdata(infiles, dataset, raw_data)
    
    # Get the raw data for the two datasets of raw data
    print('SL')
    raw_data = get_sm_rawdata(glob.glob(args.sl_infiles + '/*.xls'), 'SL', raw_data)
    print('NL')
    raw_data = get_sm_rawdata(glob.glob(args.nl_infiles + '/*.xls'), 'NL', raw_data, offset=raw_data.trap_ID.max())
    
    # Save the raw data to .csv format
    raw_data.sort_values(['dataset', 'trap_ID', 'trace']).reset_index(drop=True).to_csv(args.save_folder + '/raw_data.csv', index=False)
    minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['dataset', 'trap_ID', 'trace']).reset_index(drop=True).to_csv(args.save_folder + '/raw_data_tc.csv', index=False)
    
    # This is the order of the cycle variables in the processed dataframe
    order = phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'lineage_ID', 'generation']
    
    # The dataframe for our variables
    cycle_variables = pd.DataFrame(columns=order)
    
    for lineage_id in raw_data.lineage_ID.unique():
        print('Lineage ID:', lineage_id)
        
        # Get the lineage
        raw_lineage = raw_data[raw_data['lineage_ID'] == lineage_id]
        
        # Figure out the indices for the division events
        start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # # Check division times
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'])
        # for start, end in zip(start_indices, end_indices):
        #     plt.axvline(start, color='green')
        #     plt.axvline(end, color='red')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # the inter-division times
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        # Each cycle must consist of at least four data points
        at_least_number = np.where(cycle_durations > (2 * 0.05))[0]
        
        # Apply the four data points condition
        start_indices, end_indices, cycle_durations = start_indices[at_least_number], end_indices[at_least_number], cycle_durations[at_least_number]
        
        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / .05) + np.ones_like(cycle_durations), dtype=int)
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, int(lineage_id), fit_the_lengths=True)
        
        # Add the SM categorical variables
        cycle_variables_lineage['trap_ID'] = raw_lineage['trap_ID'].unique()[0]
        cycle_variables_lineage['trace'] = raw_lineage['trace'].unique()[0]
        cycle_variables_lineage['dataset'] = raw_lineage['dataset'].unique()[0]
        
        # Append the cycle variables to the processed dataframe
        cycle_variables = cycle_variables.append(cycle_variables_lineage[order], ignore_index=True)
    
    print('processed data:\n', cycle_variables)
    print('cleaned raw data:\n', raw_data)
    
    # reset the index for good practice
    cycle_variables.reset_index(drop=True).sort_values(['dataset', 'trap_ID', 'trace', 'generation']).to_csv(args.save_folder + '/physical_units.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['dataset', 'trap_ID', 'trace', 'generation']).to_csv(
        args.save_folder + '/trace_centered.csv', index=False)


""" Create the csv files for physical, trace-centered, and trap-centered units for MM and SM data """


def main():
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in mm_data_names:
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        args = parser.parse_args()
        
        create_folder(args.raw_data)
        create_folder(args.save_folder)
        
        main_mm(args)
        
        print('*' * 200)
    
    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time
    
    # Now we move onto the Sister Machine data
    data_origin = 'SM'
    
    parser = argparse.ArgumentParser(description='Process Sister Machine Lineage Data.')
    parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin)
    parser.add_argument('-SL', '--sl_infiles', metavar='', type=str,
                        help='Location of Sister Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/SL')
    parser.add_argument('-NL', '--nl_infiles', metavar='', type=str,
                        help='Location of Neighboring Lineage Raw Data', required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/NL')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    args = parser.parse_args()
    
    create_folder(args.raw_data)
    create_folder(args.save_folder)
    
    main_sm(args)
    
    print('*' * 200)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))


if __name__ == '__main__':
    main()
