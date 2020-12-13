#!/usr/bin/env bash

import pandas as pd
import numpy as np
import scipy.io as io
import glob
import os
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder


""" This function takes care of which statistic we want to subtract the trajectories to get what we asked for """


def minusing(info, parameters):
    tc = info.copy()
    for trap_id in info['trap_ID'].unique():
        trace = info[(info['trap_ID'] == trap_id)].copy()
        tc.loc[trace.index, parameters] = trace[parameters] - trace[parameters].mean()
                
    return tc


""" Recognizes where division events occur """


def get_division_indices(raw_trace):
    # From the raw data, we see where the difference between two values of length
    # falls drastically, suggesting a division has occurred.
    diffs = np.diff(raw_trace)
    
    # How much of a difference there has to be between two points to consider division having taken place.
    threshold_of_difference_for_division = -1.2
    
    # np.seterr(all='raise')
    
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
    
    assert len(start_indices) == len(end_indices)
    
    if raw_trace[0] in start_indices:
        del start_indices[np.where(start_indices != raw_trace[0])[0]]
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
    
    return [start_indices, end_indices]


"""We use this for linear regression in the cycle parameter process """


def linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, cycle_variables, trap_id, fit_the_lengths):
    
    # Make sure everything is aligned
    assert len(start_indices) == len(cycle_durations) == len(data_points_per_cycle) == len(end_indices)
    
    for start, end, inter_div_time, gen in zip(start_indices, end_indices, cycle_durations, np.arange(len(cycle_durations), dtype=int)):
        
        # try:
        #     assert np.rint((end-start) / 0.05) == num_of_points
        # except:
        #     print(np.rint((end-start) / 0.05))
        #     print(num_of_points)
        
        # the phenotypic variables of a cycle
        cycle = pd.Series()
        
        # For how long it grew for
        cycle['generationtime'] = inter_div_time
        
        # if len(raw_lineage['length'].iloc[start:end + 1].values) != num_of_points:
        #     print(num_of_points)
        #     print(len(raw_lineage['length']), end - start + 1)
        #     print(len(raw_lineage['length'].iloc[start:end + 1].values))
        #     exit()
        
        # for our regression
        domain = (raw_lineage['time'].iloc[start: end + 1].copy().values - raw_lineage['time'].iloc[start]).reshape(-1, 1)

        # domain = np.linspace(0, (end - start) * .05, num=(end - start + 1)).reshape(-1, 1)
        # print(domain)
        # print(np.log(raw_lineage['length'].iloc[start:end + 1].values).reshape(-1, 1))
        # print(domain)
        # print(raw_lineage['time'].iloc[start:end + 1].values)
        # exit()
        
        # domain = np.linspace(raw_lineage['time'].iloc[start], raw_lineage['time'].iloc[end], num=end - start + 1).reshape(-1, 1)
        range = np.log(raw_lineage['length'].iloc[start:end + 1].values).reshape(-1, 1)  # the end+1 is due to indexing
        
        assert len(domain) == len(range)
        
        # do the regression
        if fit_the_lengths:
            # our regression
            reg = LinearRegression().fit(domain, range)
            
            # phenotypic variables
            cycle['growth_rate'] = reg.coef_[0][0]
            cycle['length_birth'] = np.exp(reg.predict(domain[0].reshape(-1, 1))[0][0])
            cycle['length_final'] = np.exp(reg.predict(domain[-1].reshape(-1, 1))[0][0])
            cycle['fold_growth'] = cycle['growth_rate'] * cycle['generationtime']
            if start == 0:
                cycle['division_ratio'] = np.nan
            else:
                cycle['division_ratio'] = cycle['length_birth'] / cycle_variables['length_final'].iloc[-1]   # raw_lineage['length'].iloc[start] / raw_lineage['length'].iloc[start - 1]
            cycle['added_length'] = cycle['length_final'] - cycle['length_birth']
            cycle['generation'] = gen
        
        else:
            # our regression
            reg = LinearRegression().fit(domain, range)
            
            # phenotypic variables
            cycle['growth_rate'] = reg.coef_[0][0]
            cycle['length_birth'] = np.exp(range[0][0])
            cycle['length_final'] = np.exp(range[-1][0])
            cycle['fold_growth'] = cycle['growth_rate'] * cycle['generationtime']
            if start == 0:
                cycle['division_ratio'] = np.nan
            else:
                cycle['division_ratio'] = cycle['length_birth'] / cycle_variables['length_final'].iloc[-1]  # raw_lineage['length'].iloc[start] / raw_lineage['length'].iloc[start - 1]
            cycle['added_length'] = cycle['length_final'] - cycle['length_birth']
            cycle['generation'] = gen
            
        # add what lineage it is from
        cycle['trap_ID'] = trap_id
            
        # append the cycle variables to the
        cycle_variables = cycle_variables.append(cycle, ignore_index=True)
        
        # # Check the regression
        # plt.plot(domain, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.plot(domain, np.exp(range))
        # plt.show()
        # plt.close()
        
    return cycle_variables


""" Create the csv files for physical, trace-centered, and trap-centered units """


def main(args):
    print('type of MM data we are processing:', args.data_origin)
    # Directory of the MM data
    # infiles_lambda_LB = glob.glob(r'/Users/alestawsky/Downloads/MotherMachine/lambda_LB/fits/*.mat')
    # infiles_Maryam_LongTrace = glob.glob(r'/Users/alestawsky/Downloads/MotherMachine/MG1655-inLB-LongTraces/fits/*.mat')
    
    # infiles_MG1655_inLB_LongTraces = glob.glob(r'/Users/alestawsky/Downloads/MotherMachine/MG1655-inLB-LongTraces/traces/*.txt')
    # infiles_MG1655_inLB_LongTraces = glob.glob(r'/Users/alestawsky/PycharmProjects/Thesis/RawData/MG1655-inLB-LongTraces/*.txt')
    
    if args.data_origin == 'MG1655_inLB_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.txt')
    elif args.data_origin == 'Maryam_LongTraces':
        infiles = glob.glob(args.raw_data + '/*.csv')
        infiles = infiles + glob.glob(args.raw_data + '/*.xls')
    else:
        infiles = glob.glob(args.raw_data + '/*.txt')
        
    # print('infiles_MG1655_inLB_LongTraces', infiles, sep='\n')

    # the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['trap_ID', 'generation'])
    
    # The dataframe for our raw data lineages
    raw_data = pd.DataFrame(columns=['time', 'length', 'trap_ID'])
    
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
        '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos17-1.txt',
        
    ]

    # load first sheet of each Excel-File, fill internal data structure]
    for count, filename in enumerate(infiles):
        
        print(count, filename.split('/')[-1], sep=': ')
        
        # if filename in filenames_with_nans:
        #     continue
    
        # creates a dataframe from the .txt or .csv file
        if args.data_origin == 'MG1655_inLB_LongTraces':
            if filename.split('/')[-1] == 'pos4-4.txt':
                print('In this particular case the lineage divides after the first time-point and it has an extra column.')
                # This is because in this particular case the lineage divides after the first time-point and it has an extra column
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']].reset_index(
                    drop=True)
                raw_lineage = raw_lineage.iloc[1:].reset_index(drop=True)
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args.data_origin == 'Maryam_LongTraces':
            try:
                raw_lineage = pd.read_csv(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
            except:
                raw_lineage = pd.read_excel(filename, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
        elif args.data_origin == 'lambda_LB':
            if filename in extra_column:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            else:
                raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        elif args.data_origin == 'LAC_M9':
            raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        
        # Make the time-steps accurate to two decimal points
        raw_lineage['time'] = raw_lineage['time'].round(2)
        
        # the trap ID
        raw_lineage['trap_ID'] = count
        
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        if not np.all(step_sizes == step_sizes[0]):
            print(filename, ' has steps that are not the same size, are we are not gonna use these...')
            continue
        
        # Make sure there are no NaNs
        if raw_lineage.isna().values.any():
            print(raw_lineage.isna().values.any())
            print(raw_lineage.isna().sum())
            print(raw_lineage)
            print(raw_lineage[raw_lineage.isnull()])
            exit()
        
        # if filename == '/Users/alestawsky/PycharmProjects/Thesis/RawData/lambda_LB/pos17-1.txt':
        #     print(raw_lineage)
        #     exit()
        
        # print(raw_lineage)
    
        # add it to the total data
        raw_data = raw_data.append(raw_lineage, ignore_index=True)
        
        # Figure out the indices for the division events
        start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        # if (args.data_origin == 'MG1655_inLB_LongTraces') & (filename == args.raw_data + '/pos4-4.txt'):
        #     # This is because in this particular case the lineage divides after the first time-point
        #     # start_indices = start_indices[1:]
        #
        #     # Check division times
        #     plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'])
        #     for start, end in zip(start_indices, end_indices):
        #         plt.axvline(start, color='green')
        #         plt.axvline(end, color='red')
        #     plt.tight_layout()
        #     plt.show()
        #     plt.close()

        assert len(start_indices) == len(end_indices)

        # the inter-division times
        cycle_durations = raw_lineage['time'].values[end_indices] - raw_lineage['time'].values[start_indices]
        
        # Each cycle must consist of at least four data points
        at_least_number = np.where(cycle_durations > (2 * 0.05))[0]
        
        # print('condition', len(cycle_durations), len(cycle_durations[at_least_number]))
        
        start_indices, end_indices, cycle_durations = start_indices[at_least_number], end_indices[at_least_number], cycle_durations[at_least_number]
        
        assert len(start_indices) == len(end_indices) == len(cycle_durations)

        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / .05) + np.ones_like(cycle_durations), dtype=int)

        # add the cycle variables to the overall dataframe
        cycle_variables = linear_regression(raw_lineage, cycle_durations, start_indices, end_indices, data_points_per_cycle, cycle_variables, int(count), fit_the_lengths=True)
        
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
    raw_data.reset_index(drop=True).to_csv(args.save_folder+'/raw_data.csv')
    cycle_variables.reset_index(drop=True).to_csv(args.save_folder+'/physical_units.csv')
    minusing(raw_data.reset_index(drop=True), ['length']).to_csv(args.save_folder+'/raw_data_tc.csv')
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).to_csv(args.save_folder+'/trace_centered.csv')


if __name__ == '__main__':
    import argparse
    import os
    import time

    # How long does running this take?
    first_time = time.time()
    
    for data_origin in ['lambda_LB', 'MG1655_inLB_LongTraces', 'Maryam_LongTraces', 'lambda_LB', 'LAC_M9']:
        
        parser = argparse.ArgumentParser(description='Process Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-raw_data', '--raw_data', metavar='', type=str, help='Raw Data location.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/RawData/' + data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                            required=False, default='trace_centered.csv')
        parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        args = parser.parse_args()
        
        create_folder(args.raw_data)
        create_folder(args.save_folder)
        
        main(args)

        print('*' * 200)

    print("--- took %s minutes ---" % ((time.time() - first_time) / 60))
