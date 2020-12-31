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
    
    # make them integers
    cycle_variables_lineage['lineage_ID'] = int(lin_id)
    cycle_variables_lineage['generation'] = np.arange(len(cycle_variables_lineage), dtype=int)
    cycle_variables_lineage = cycle_variables_lineage.sort_values('generation')
    
    return [cycle_variables_lineage.reset_index(drop=True), without_nans]


""" Create the csv files for physical, trace-centered, and trap-centered units for MM data of a certain type """


def start_process(args):
    print('type of MM data we are processing:', args['data_origin'])
    
    infiles = glob.glob(args['raw_data'] + '/*')
    # filenames = [file.split('/')[-1].split('.')[0] for file in infiles]
    # extensions = [file.split('/')[-1].split('.')[1] for file in infiles]
    # print(infiles)
    # print(filenames)
    # print(extensions)
    # print(len(infiles), len(filenames), len(extensions))
    # exit()
    
    # the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    with_outliers_cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    # The dataframe for our raw data lineages
    raw_data = pd.DataFrame(columns=['time', 'length', 'lineage_ID', 'filename'])
    
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
    
    # When we take out the outliers we keep trace of how many cycles we lost wrt. the lineage and the pooled ensemble
    whats_left = {'variable': [], 'lineage': [], 'number_taken_out': []}
    
    # In case we can't use some files we want the lineage IDs to be in integer order
    offset = 0
    
    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(infiles):
    
        filename = file.split('/')[-1].split('.')[0]
        extension = file.split('/')[-1].split('.')[1]
        
        # Tells us the trap ID and the source (filename)
        print(count, filename, sep=': ')
        
        if extension == 'txt':
            if filename.split('/')[-1].split('.')[0] in extra_column:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            elif filename.split('/')[-1] == 'pos15.txt':
                print('This one is special.')
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein', '__', '___', '___1'])[
                    ['time', 'length']]
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif filename.split('/')[-1].split('.')[1] == 'csv':
            pass
        elif filename.split('/')[-1].split('.')[1] == 'xls':
            pass
        else:
            raise IOError('Data extension not supported')
        
        if args.data_origin == 'Maryam_LongTraces':
            # This is because some the data is in .xls format while others are in .csv
            try:
                raw_lineage = pd.read_csv(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
            except:
                raw_lineage = pd.read_excel(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        elif args.data_origin == 'lambda_LB':
            # There are quite a lot of files with an extra column at the beginning
            if filename in extra_column:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            elif filename.split('/')[-1] == 'pos15.txt':
                print('This one is special.')
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein', '__', '___', '___1'])[
                    ['time', 'length']]
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time'].round(2)
        raw_lineage['filename'] = filename.split('/')[-1]
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
        cycle_variables_lineage, with_outliers = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, int(count - offset), fit_the_lengths=True)
        
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
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data.csv', index=False)
    minusing(raw_data.reset_index(drop=True), ['length']).sort_values(['lineage_ID']).to_csv(args.save_folder + '/raw_data_tc.csv', index=False)
    cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args.save_folder + '/physical_units.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args.save_folder + '/trace_centered.csv',
                                                                                                                                                   index=False)
    
    with_outliers_cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args.save_folder + '/physical_units_with_outliers.csv', index=False)
    minusing(with_outliers_cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        args.save_folder + '/trace_centered_with_outliers.csv',
        index=False)


""" Create the csv files for physical, trace-centered, and trap-centered units for MM and SM data """


def main():
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in dataset_names + ['SM']:
        print(data_origin)
        
        args = {
            'data_origin': data_origin,
            'raw_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin,
            'save_folder': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin,
            # The following two are only needed for SM data (Sister and Neighbor Lineages) but we put it in anyways for cleanliness #
            'sl': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/SL',
            'nl': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/SM/NL'
        }
        
        # Make sure the folders where we place the data are created already
        create_folder(args['raw_data'])
        create_folder(args['save_folder'])
        
        # compare_cycle_variables_to_raw_data(args)
        
        start_process(args)
        
        print('*' * 200)
        print('took {:.2}'.format((time.time() - first_time)/60)+' minutes')
        print('*' * 200)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))


if __name__ == '__main__':
    main()
