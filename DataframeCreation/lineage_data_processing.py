#!/usr/bin/env bash

# import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.linear_model import LinearRegression
# import argparse
from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder
# sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
# from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder


"""We use this for linear regression in the cycle parameter process """


def linear_regression(data, dataID, discretize_by, ks, cycle_durations, start_indices, end_indices, data_points_per_cycle, fit_the_length_at_birth):
    # The x and y values per generation for the regression to get the growth rate
    domains_for_regression_per_gen = [np.linspace(data[dataID]['time' + ks][start], data[dataID]['time' + ks][end], num=num_ofdata_points) for start, end, num_ofdata_points in
                                      zip(start_indices, end_indices, data_points_per_cycle)]
    ranges_for_regression_per_gen = [np.log(data[dataID][discretize_by + ks][start:end + 1])  # the end+1 is due to indexing
                                     for start, end in zip(start_indices, end_indices)]
    
    # Arrays where the intercepts (the length at birth) and the slopes (the growth rate) will be stored
    regression_intercepts_per_gen = []
    regression_slopes_per_gen = []
    for domain, y_vals in zip(domains_for_regression_per_gen, ranges_for_regression_per_gen):  # loop over generations in the trace
        # reshape the x and y values
        domain = np.array(domain).reshape(-1, 1)
        y_vals = np.array(y_vals).reshape(-1, 1)
        # do the regression
        reg = LinearRegression().fit(domain, y_vals)
        # save them to their respective arrays
        regression_slopes_per_gen.append(reg.coef_[0][0])
        regression_intercepts_per_gen.append(np.exp(reg.predict(domain[0].reshape(-1, 1))[0][0]))
    
    # change them to numpy arrays
    regression_slopes_per_gen = np.array(regression_slopes_per_gen)
    regression_intercepts_per_gen = np.array(regression_intercepts_per_gen)
    
    # check if the growth_rate is negative, meaning it is obviously a not a generation
    checking_pos = np.array([slope > 0 for slope in regression_slopes_per_gen])
    if not checking_pos.all():
        print('----')
        print("there's been a negative or zero growth_rate found!")
        print('{} {}'.format(ks, dataID))
        print('----')
        # change growth_rate
        regression_slopes_per_gen = regression_slopes_per_gen[np.where(checking_pos)]
        # change length at birth
        regression_intercepts_per_gen = regression_intercepts_per_gen[np.where(checking_pos)]
        # change generationtime
        cycle_durations = cycle_durations[np.where(checking_pos)]
        # change the indices which will change the rest
        start_indices, end_indices = start_indices[np.where(checking_pos)], end_indices[np.where(checking_pos)]
    
    # check if the length_final <= length_birth, meaning it is obviously a not a generation
    if fit_the_length_at_birth:
        checking_pos = np.array(
            [length_birth < length_final for length_birth, length_final in zip(regression_intercepts_per_gen, np.array(data[dataID][discretize_by + ks][end_indices]))])
    else:
        checking_pos = np.array([length_birth < length_final for length_birth, length_final in
                                 zip(np.array(data[dataID][discretize_by + ks][start_indices]), np.array(data[dataID][discretize_by + ks][end_indices]))])
    
    # If we have any of these problems then we disregard the generation
    if not checking_pos.all():
        print('length_final <= length_birth, and this so-called "generation" was taken out of its dataframe')
        # change growth_rate
        regression_slopes_per_gen = regression_slopes_per_gen[np.where(checking_pos)]
        # change length at birth
        regression_intercepts_per_gen = regression_intercepts_per_gen[np.where(checking_pos)]
        # change generationtime
        cycle_durations = cycle_durations[np.where(checking_pos)]
        # change the indices which will change the rest
        start_indices, end_indices = start_indices[np.where(checking_pos)], end_indices[np.where(checking_pos)]
    
    return [cycle_durations, start_indices, end_indices, regression_slopes_per_gen, regression_intercepts_per_gen]


""" This function takes care of which statistic we want to subtract the trajectories to get what we asked for """


def minusing(Traj_A, Traj_B, kind, parameters):
    if kind != 'physical_units':
        # we want to subtract the time-average to each value
        if kind == 'trap_centered':
            trap = np.sum(pd.concat([Traj_A, Traj_B], axis=0)[parameters]) / len(pd.concat([Traj_A, Traj_B], axis=0)[parameters])
            Traj_A[parameters] = Traj_A[parameters] - trap[parameters]
            Traj_B[parameters] = Traj_B[parameters] - trap[parameters]
        elif kind == 'trace_centered':
            Traj_A[parameters] = Traj_A[parameters] - np.sum(Traj_A[parameters]) / len(Traj_A[parameters])
            Traj_B[parameters] = Traj_B[parameters] - np.sum(Traj_B[parameters]) / len(Traj_B[parameters])
    
    return Traj_A, Traj_B


""" check division times, ie. plot all the A and B trajectories with the division times highlighted """


def division_times_check(**kwargs):
    # ask if we want to save or show the division check points
    answer = input("Do you want to save the plots or do you want to see them?")
    if answer == 'save':
        name = input("What is the name of the file they should be put under?")
        create_folder(name)
    elif answer == 'see':
        pass
    else:
        print('wrong input! you will only see them.')
    
    # The observable we will use to see when the cell divided and what the cycle data is based on
    discretize_by = kwargs.get('discretize_by', [])
    
    # which initial and final generations to include in the data
    start_index = kwargs.get('start_index', None)
    end_index = kwargs.get('end_index', None)
    dataset = kwargs.get('dataset', '')
    
    # where to put the raw data
    raw_data = list()
    keylist = list()
    
    # Get the keylists for sisters pairs and neighbor pairs whos' raw data are seperate in the folders attached
    raw_data, keylist = get_keylist(raw_data, keylist, infiles=kwargs.get('infiles', []))
    
    for raw_dataID in range(len(raw_data)):
        for ks in ['A', 'B']:
            # plot the trajectory and the beginning and end points
            plt.plot(raw_data[raw_dataID]['time' + ks], raw_data[raw_dataID][discretize_by + ks], marker='+')
            start_indices, end_indices = get_division_indices(data=raw_data, data_id=raw_dataID, discretize_by=discretize_by, ks=ks)
            plt.scatter(raw_data[raw_dataID]['time' + ks].iloc[start_indices], raw_data[raw_dataID][discretize_by + ks].iloc[start_indices], color='red', label='first point in cycle')
            plt.scatter(raw_data[raw_dataID]['time' + ks].iloc[end_indices], raw_data[raw_dataID][discretize_by + ks].iloc[end_indices], color='green', label='last point in cycle')
            plt.legend()
            
            # save them or show them based on what we answered
            if answer == 'save':
                plt.savefig(name + '/' + ks + '_' + str(raw_dataID) + '.png', dpi=300)
            elif answer == 'see':
                plt.show()
            else:
                plt.show()
            plt.close()


""" Use this to get the starting and ending indices in the raw data """


def get_division_indices(data, data_id, discretize_by, ks):
    # From the raw data, we see where the difference between two values of 'discretize_by' falls drastically,
    # suggesting a division has occurred.
    diffdata = np.diff(data[data_id][discretize_by + ks])
    
    # How much of a difference there has to be between two points to consider division having taken place.
    threshold_of_difference_for_division = -1
    
    np.seterr(all='raise')
    
    # The array of indices where division took place.
    index_div = np.where(diffdata < threshold_of_difference_for_division)[0].flatten()
    
    # If the two consecutive indices in the array of indices where division took place are less than two time steps away, then we discard
    # them
    for ind1, ind2 in zip(index_div[:-1], index_div[1:]):
        if ind2 - ind1 <= 2:
            index_div = np.delete(index_div, np.where(index_div == ind2))
    
    # An array of indices where a cycle starts and ends
    start_indices = [x + 1 for x in index_div]
    end_indices = [x for x in index_div]
    start_indices.append(0)  # to count the first cycle
    start_indices.sort()  # because '0' above was added at the end and the rest are in order already
    del start_indices[-1]  # we do not count the last cycle because it most likely did not finish by the time recording ended
    end_indices.sort()  # same as with start_indices, just in case
    
    # If the last index is in the array of indices where a cycle starts, then it obviously a mistake and must be removed
    if data[data_id][discretize_by + ks].index[-1] in start_indices:
        start_indices.remove(data[data_id][discretize_by + ks].index[-1])
    
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
    
    return start_indices, end_indices


"""
Transform measured time series data into data for each generation:
(1) Find cell division events as a large enough drop in the measured value given by 'discretize_by'
(2) Compute various observables from these generations of cells
(3) Returns two pandas dataframes for each of the discretized trajectories
"""


def cycle_parameter_process(**kwargs):
    # arbitrary but necessary specification options for the two cells inside the trap
    keysuffix = ['A', 'B']
    
    dataID, data, discretize_by, fit_the_length_at_birth, dataset = kwargs['dataID'], kwargs['data'], kwargs['discretize_by'], kwargs['fit_the_length_at_birth'], kwargs['dataset']
    
    # Return two pandas-dataframes for the two trajectories usually contained in one file as two elements of a list.
    # As the two sisters do not need to have the same number of cell divisions,
    # a single dataframe might cause problems with variably lengthed trajectories.
    ret = list()
    
    for ks in keysuffix:
        
        # use this to get the starting and ending indices in the raw data
        start_indices, end_indices = get_division_indices(data_id=dataID, data=data, discretize_by=discretize_by, ks=ks)
        
        # How long each cycle is
        cycle_durations = np.array(data[dataID]['time' + ks][end_indices]) - np.array(data[dataID]['time' + ks][start_indices])
        
        # Store results in this dictionary, which can be easier transformed into pandas-dataframe
        ret_ks = dict()
        
        # Number of raw data points per generation/cycle
        data_points_per_cycle = np.array(np.rint(cycle_durations / .05) + np.ones_like(cycle_durations), dtype=int)
        
        # Do the linear regression and check if it is really a generation, ie. length_final > length_birth, if growth_rate > 0
        cycle_durations, start_indices, end_indices, regression_slopes_per_gen, regression_intercepts_per_gen = linear_regression(data, dataID, discretize_by, ks, cycle_durations, start_indices,
                                                                                                                                  end_indices,
                                                                                                                                  data_points_per_cycle, fit_the_length_at_birth)
        
        # Duration of growth before division
        ret_ks['generationtime'] = np.around(cycle_durations, decimals=2)
        
        # Due to limitations of data, mainly that after some obvious division points the length of the bacteria drops, which means that it is
        # shrinking, which we assume is impossible. Either way we present is as an option, just in case.
        if fit_the_length_at_birth:
            ret_ks['length_birth'] = regression_intercepts_per_gen
        else:
            ret_ks['length_birth'] = np.array(data[dataID][discretize_by + ks][start_indices])
        
        # The measured final length before division, however it is worth to note that we have a very good approximation of this observable by
        # the mapping presented in Susman et al. (2018)
        ret_ks['length_final'] = np.array(data[dataID][discretize_by + ks][end_indices])
        
        # The rate at which the bacteria grows in this cycle. NOTICE THAT THIS CHANGED FROM 'growth_length'
        ret_ks['growth_rate'] = regression_slopes_per_gen
        
        # The fold growth, ie. how much it grew. Defined by the rate at which it grew multiplied by the time it grew for
        # NOTE: We call this 'phi' before, in an earlier version of the code
        ret_ks['fold_growth'] = ret_ks['generationtime'] * ret_ks['growth_rate']
        
        # Calculating the division ratios, percentage of mother length that went to daughter
        div_rats = []
        for final, next_beg in zip(ret_ks['length_final'][:-1], ret_ks['length_birth'][1:]):
            div_rats.append(next_beg / final)
        
        # we use the length at birth that is not in the dataframe in order to get enough division ratios
        div_rats.append(data[dataID][discretize_by + ks][end_indices[-1] + 1] / ret_ks['length_final'][-1])
        ret_ks['division_ratio'] = div_rats
        
        # the added length to check the adder model
        ret_ks['added_length'] = ret_ks['length_final'] - ret_ks['length_birth']
        
        # Make into dataframe
        ret_ks = pd.DataFrame(ret_ks)
        
        # Categorical labels
        ret_ks['dataset'] = dataset
        ret_ks['trap_ID'] = str(dataID)
        ret_ks['trace'] = ks
        ret_ks['generation'] = np.arange(len(ret_ks))
        
        # we have everything, now make a dataframe
        ret.append(ret_ks)
    return ret


""" Reading the raw data into an array of pd.Dataframes """


def get_keylist(data, keylist, infiles):
    # load first sheet of each Excel-File, fill internal data structure]
    for filename in infiles:
        
        # creates a dataframe from the excel file
        tmpdata = pd.read_excel(filename).dropna(axis='index', how='any')
        
        # This is because Neighbor Cells raw csv files contain the columns L1, L2 in place for lengthA, lengthB sometimes
        if 'lengthA' not in tmpdata.columns:
            # creates a dataframe from the excel file
            tmpdata = pd.read_excel(filename).dropna(axis='index', how='any').rename(columns={'L1': 'lengthA', 'L2': 'lengthB'})
        
        # the data contains all dataframes from the excel files in the directory _infiles
        data.append(tmpdata)
        for k in tmpdata.keys():
            if not str(k) in keylist:
                # this list contains the column names of the excel file
                keylist.append(str(k))
    
    # There must be some data, just a check
    assert len(data) > 0
    
    return [data, keylist]


""" adds each trace to the big information dataframe """


def create_information_dataframe(**kwargs):
    # The observable we will use to see when the cell divided and what the cycle data is based on
    discretize_by = kwargs.get('discretize_by', [])
    
    # which initial and final generations to include in the data
    start_index = kwargs.get('start_index', None)
    end_index = kwargs.get('end_index', None)
    dataset = kwargs.get('dataset', '')
    
    # where to put the raw data
    raw_data = list()
    keylist = list()
    
    # Get the keylists for sisters pairs and neighbor pairs whos' raw data are seperate in the folders attached
    raw_data, keylist = get_keylist(raw_data, keylist, infiles=kwargs.get('infiles', []))
    
    # info = pd.DataFrame(columns=phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'generation'])
    info = kwargs.get('info_df')
    
    # loop over all traces
    for experiment_number in range(len(raw_data)):
        traj_a, traj_b = cycle_parameter_process(dataID=experiment_number,
                                                 discretize_by=discretize_by,
                                                 fit_the_length_at_birth=False,
                                                 dataset=dataset,
                                                 data_len=len(raw_data),
                                                 keylist=keylist,
                                                 data=raw_data,
                                                 start_index=start_index,
                                                 end_index=end_index)
        
        # apply the hard generation limits across all traces
        traj_a = traj_a.iloc[start_index:end_index].reset_index(drop=True)
        traj_b = traj_b.iloc[start_index:end_index].reset_index(drop=True)
        
        # If we need to subtract by trace or trap time-average we do so here
        traj_a, traj_b = minusing(traj_a, traj_b, kwargs['kind'], phenotypic_variables)
        
        # save these two trajectories to the information dataframe
        info = info.append(traj_a.append(traj_b, ignore_index=True), ignore_index=True).reset_index(drop=True)
    
    return info


def main(args):
    #sp_infiles=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/SP',
    # nc_infiles=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/NC',
    # save_folder=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data'
    
    # visually_check_divisions()
    
    print('-' * 50)
    print('-' * 50)
    
    # For physical_units and trace_centered unless any one otherwise specified
    for k, name in zip(['physical_units', 'trace_centered'], [args.pu, args.tc]):
        
        info_df = pd.DataFrame(columns=phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'generation'])
        
        # For every type of lineage we want to look at. These are all the options SP and NC, not CTRL, hence [:-1]
        for dataset, infiles in zip(datasets[:-1], [glob.glob(args.sp_infiles + '/*.xls'), glob.glob(args.nc_infiles + '/*.xls')]):
            info_df = create_information_dataframe(infiles=infiles, discretize_by='length', kind=k, start_index=None, end_index=None, dataset=dataset, info_df=info_df)
        
        # Save the dataframe to a csv
        info_df.to_csv('{}/{}'.format(args.save_folder, name), index=False)
        
        print('-' * 50)
        print('-' * 50)

# # visually_check_divisions()
# parser = argparse.ArgumentParser(description='Process Lineage Data.')
# parser.add_argument('-SP', '--sp_infiles', metavar='', type=str,
#                     help='Location of Sister Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/SP')
# parser.add_argument('-NC', '--nc_infiles', metavar='', type=str,
#                     help='Location of Neighboring Cell Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/NC')
# parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
# # group = parser.add_mutually_exclusive_group()
# # group.add_argument('-de', '--division_events', action='store_true', help='0: Do not show the division events.\n 1: Show the division events.')
#
# args = parser.parse_args()
#
# create_folder(args.save_folder)
#
# print('-' * 50)
# print('-' * 50)
#
# # For physical_units and trace_centered unless any one otherwise specified
# for k in ['physical_units', 'trace_centered']:
#     print(k)
#     print('-' * 50)
#     print('-' * 50)
#
#     info_df = pd.DataFrame(columns=phenotypic_variables + ['dataset', 'trap_ID', 'trace', 'generation'])
#
#     # For every type of lineage we want to look at. These are all the options SP and NC, not CTRL, hence [:-1]
#     for dataset, infiles in zip(datasets[:-1], [glob.glob(args.sp_infiles + '/*.xls'), glob.glob(args.nc_infiles + '/*.xls')]):
#         print(dataset)
#         info_df = create_information_dataframe(infiles=infiles, discretize_by='length', kind=k, start_index=None, end_index=None, dataset=dataset, info_df=info_df)
#
#     # Save the dataframe to a csv
#     info_df.to_csv('{}/{}.csv'.format(args.save_folder, k), index=False)
#
#     print('-' * 50)
#     print('-' * 50)


if __name__ == '__main__':
    import os
    import argparse
    import time
    
    print('running from __main__')

    first_time = time.time()
    
    parser = argparse.ArgumentParser(description='Process Lineage Data.')
    parser.add_argument('-SP', '--sp_infiles', metavar='', type=str,
                        help='Location of Sister Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/SP')
    parser.add_argument('-NC', '--nc_infiles', metavar='', type=str,
                        help='Location of Neighboring Cell Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/NC')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Data')
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                        required=False, default='trace_centered.csv')

    args = parser.parse_args()
    
    main(args)

    print("--- %s seconds ---" % (time.time() - first_time))
