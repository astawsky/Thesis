#!/usr/bin/env bash

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from scipy.stats import zscore
from sklearn.mixture import GaussianMixture
from scipy.signal import find_peaks, peak_prominences
from sklearn.linear_model import LinearRegression
from AnalysisCode.global_variables import phenotypic_variables, create_folder, dataset_names, sm_datasets, mm_datasets, wang_datasets, tanouchi_datasets, seaborn_preamble, cgsc_6300_wang_exps

""" This function takes care of which statistic we want to subtract the trajectories to get what we asked for """


def minusing(info, parameters):
    tc = info.copy()
    for lineage_id in info['lineage_ID'].unique():
        trace = info[(info['lineage_ID'] == lineage_id)].copy()
        tc.loc[trace.index, parameters] = trace[parameters] - trace[parameters].mean()
    
    return tc


""" Check if there are any "singular" points due to camera problems maybe and eliminate them """


def eliminate_drop_offs(raw_lineage):
    np.diff(raw_lineage)


""" Recognizes where division events occur """


def get_division_indices(raw_trace):
    # From the raw data, we see where the difference between two values of length
    # falls drastically, suggesting a division has occurred.
    diffs = -np.diff(np.log(raw_trace[np.where(~np.isnan(raw_trace))]))
    
    peaks, _ = find_peaks(diffs, threshold=np.log(1.3))
    
    start_indices = np.append([0], peaks[:-1] + 1)
    end_indices = peaks
    
    # If the first cycle is too small to be a cycle
    if start_indices[1] - start_indices[0] < 5:
        start_indices = start_indices[1:]
        end_indices = end_indices[1:]
    
    # Make sure they are the same size
    assert len(start_indices) == len(end_indices)
    
    # plt.plot(raw_trace[np.where(~np.isnan(raw_trace))])
    # plt.scatter(start_indices, raw_trace[np.where(~np.isnan(raw_trace))][start_indices], color='green')
    # plt.scatter(end_indices, raw_trace[np.where(~np.isnan(raw_trace))][end_indices], color='red')
    # # plt.axvline()
    # plt.show()
    # plt.close()
    
    return [start_indices, end_indices]


"""We use this for linear regression in the cycle parameter process """


def linear_regression(clean_lin, raw_lineage, start_indices, end_indices, lin_id, fit_the_lengths):
    # the dataframe for our variables
    cycle_variables_lineage = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    if clean_lin.isnull().values.any():
        print('NaNs found!')
        print(clean_lin[clean_lin.isnull()])
        exit()
    
    rl = clean_lin[['time', 'length']].copy().dropna()
    
    assert len(start_indices) == len(end_indices)
    
    new_starting = []
    new_ending = []
    for start, end, gen in zip(start_indices, end_indices, np.arange(len(start_indices), dtype=int)):
        
        # for our regression
        # domain1 = np.linspace(np.min(rl['time'].iloc[start: end].copy().values - rl['time'].iloc[start]),
        #                      np.max(rl['time'].iloc[start: end].copy().values - rl['time'].iloc[start]),
        #                      num=int(np.round(np.max(rl['time'].iloc[end] - rl['time'].iloc[start]) / step_size, 0))).reshape(-1, 1)
        # domain2 = np.arange(
        #     0.0, np.max(rl['time'].iloc[start: end].copy().values - rl['time'].iloc[start]), step_size
        # ).reshape(-1, 1)
        
        if end > len(rl):
            print('used the new end!')
            exit()
            end_indices = np.where(end_indices == end, len(rl) - 1, end_indices)
            end = len(rl) - 1
        domain = (clean_lin['time'].iloc[start: end + 1].copy().values - clean_lin['time'].iloc[start]).reshape(-1, 1)
        
        if len(domain) == 0:
            print('Only NaNs')
            # Take this out
            start_indices = start_indices[np.where(start_indices != start)]
            end_indices = end_indices[np.where(end_indices != end)]
            continue
        
        # print(domain)
        # print('{}'*100)
        
        # if np.max(domain) != np.max(rl['time'].iloc[start: end].copy().values - rl['time'].iloc[start]):
        #     print(np.max(domain), np.max(rl['time'].iloc[start: end].copy().values - rl['time'].iloc[start]))
        
        # if ppc != len(domain):
        #     print(ppc, len(domain))
        #     raise IOError('points per cycle and length of domain are not the same! {} != {}'.format(ppc, len(domain)))
        
        # domain = np.linspace(clean_lin['time'].iloc[start], clean_lin['time'].iloc[end], num=end - start + 1).reshape(-1, 1)
        range = np.log(rl['length'].iloc[start:end + 1].values).reshape(-1, 1)  # the end+1 is due to indexing
        
        # Make sure they are the same size for the regression!
        if len(domain) != len(range):
            print(len(clean_lin[['time', 'length']].iloc[start: end]))
            print(len(domain), len(range))
            exit()
        
        if (rl['length'].iloc[start:end].values <= 0).any():
            print(rl['length'].iloc[start:end])
            plt.plot(np.arange(start, end), rl['length'].iloc[start:end], ls='--')
            plt.plot(rl['length'])
            plt.show()
            plt.close()
        
        # our regression
        try:
            reg = LinearRegression().fit(domain, range)
        except:
            print('exception!')
            print(domain)
            print(range)
            exit()
        
        # If the growth rate is non-positive then it is obviously not a credible cycle, probably a glitch
        if (reg.coef_[0][0] <= 0):
            print('Negative Growth rate in cycle! Not counted! {}'.format(reg.coef_[0][0]))
            # Take this out
            start_indices = start_indices[np.where(start_indices != start)]
            end_indices = end_indices[np.where(end_indices != end)]
            
            # plt.plot(np.arange(start, end), rl['length'].iloc[start:end], ls='--')
            # plt.plot(rl['length'])
            # plt.show()
            # plt.close()
            continue
        if (clean_lin.iloc[start:end + 1].count()['length'] <= ((2 / 3) * len(clean_lin.iloc[start:end + 1]))):
            print('A lot of NaNs! Not counted!')
            # Take this out
            start_indices = start_indices[np.where(start_indices != start)]
            end_indices = end_indices[np.where(end_indices != end)]
            continue
        
        # the phenotypic variables of a cycle
        cycle = pd.Series()
        
        # Define the cycle phenotypic variables
        cycle['generationtime'] = rl['time'].iloc[end + 1] - rl['time'].iloc[start]
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
        
        # raw_start = np.where(raw_lineage['length'].values == clean_lin['length'].iloc[start])[0]
        # raw_end = np.where(raw_lineage['length'].values == clean_lin['length'].iloc[end])[0]
        #
        # # print(raw_start, start)
        # # print(raw_end, end)
        #
        # print('raw_start', raw_start, 'raw_end', raw_end, sep='\n')
        # # print(end_indices)
        #
        # if len(raw_start) == 1:
        #     raw_start = raw_start[0]
        # else:
        #     raw_start = raw_start[np.where(raw_start < new_starting[-1])]
        #     # raw_start = raw_start[[~any([False if (s < n < e) else True for s, e in zip(start_indices, end_indices)]) for n in raw_start]][0]
        #
        # if len(raw_end) == 1:
        #     raw_end = raw_end[0]
        # else:
        #     raw_end = raw_end[np.where(raw_end < new_ending[-1])]
        #     # raw_end = raw_end[[~any([False if (s < n < e) else True for s, e in zip(start_indices, end_indices)]) for n in raw_end]][0]
        #
        #     # newn = np.setdiff1d(raw_end, start_indices)
        #     # raw_end = newn[[~any([False if (s < n < e) else True for s, e in zip(start_indices, end_indices)]) for n in newn]]
        #
        # print('raw_start', raw_start, 'raw_end', raw_end, sep='\n')
        # print('-'*200)
        
        # raw_start = raw_start[0] if len(raw_start) == 1 else np.intersect1d(raw_start, start_indices)[0]
        # raw_end = raw_end[0] if len(raw_end) == 1 else np.intersect1d(raw_end, end_indices)[0]
        
        # raw_start = raw_start[0] if len(raw_start) == 1 else raw_start[np.intersect1d(raw_start, new_starting, return_indices=True)[1]][0]
        # raw_end = raw_end[0] if len(raw_end) > 1 else raw_end[np.intersect1d(raw_end, new_ending, return_indices=True)[1]][0]
        
        # np.setdiff1d(raw_start, new_starting)[0]
        
        # new_starting.append(raw_start)
        # new_ending.append(raw_end)
        
        # print(start)
        # print(cycle)
        # print(cycle_variables_lineage)
        
        # Check the regression
        # print(domain)
        # print(clean_lin['time'].iloc[start:end+1].values)
        # plt.plot(domain, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.plot(domain, np.exp(range))
        
        # plt.plot(clean_lin['time'].values, clean_lin['length'].values, color='orange')
        # plt.scatter(clean_lin['time'].iloc[start_indices].values, clean_lin['length'].iloc[start_indices].values, color='blue')
        # plt.plot(clean_lin['time'].iloc[start:end].values, [cycle['length_birth'] * np.exp(cycle['growth_rate'] * dom) for dom in domain])
        # plt.scatter(clean_lin['time'].iloc[start], cycle['length_birth'])
        # plt.yscale('log')
        # plt.show()
        # plt.close()
    
    # The experimental size variables
    # cycle_variables_lineage['div_then_fold'] = np.append(np.nan, cycle_variables_lineage['division_ratio'].values[:-1] * np.exp(cycle_variables_lineage['fold_growth'].values[1:]))
    cycle_variables_lineage['div_and_fold'] = cycle_variables_lineage['division_ratio'] * cycle_variables_lineage['fold_growth']
    # cycle_variables_lineage['fold_then_div'] = np.append(cycle_variables_lineage['division_ratio'].values[1:] * np.exp(cycle_variables_lineage['fold_growth'].values[:-1]), np.nan)
    
    # Good practice
    cycle_variables_lineage = cycle_variables_lineage.sort_values('generation')
    
    # Without throwing away outliers
    without_nans = cycle_variables_lineage.copy().reset_index(drop=True)
    
    insider_condition = np.abs(cycle_variables_lineage[phenotypic_variables] - cycle_variables_lineage[phenotypic_variables].mean()) < (3 * cycle_variables_lineage[phenotypic_variables].std())
    
    # Throwing away outliers
    cycle_variables_lineage[phenotypic_variables] = cycle_variables_lineage[phenotypic_variables].where(
        insider_condition,
        other=np.nan
    )
    
    # Without outliers
    outlier_start_indices = start_indices[[True if row.any() else False for row in ~insider_condition.values]]
    outlier_end_indices = end_indices[[True if row.any() else False for row in ~insider_condition.values]]  # [~insider_condition]
    
    # make them integers
    cycle_variables_lineage['lineage_ID'] = int(lin_id)
    cycle_variables_lineage['generation'] = np.arange(len(cycle_variables_lineage), dtype=int)
    cycle_variables_lineage = cycle_variables_lineage.sort_values('generation')
    
    # print(new_starting, new_ending, sep='\n')
    
    # plt.plot(raw_lineage['length'].values)
    # plt.scatter(new_starting, raw_lineage['length'].values[new_starting], label='start', c='green')
    # plt.scatter(new_ending, raw_lineage['length'].values[new_ending], label='end', c='red')
    # plt.show()
    # plt.close()
    
    return [cycle_variables_lineage.reset_index(drop=True), without_nans, start_indices, end_indices, outlier_start_indices, outlier_end_indices]  # start_indices, end_indices


""" Gets rid of the machine error from the signal as well as possible """


def clean_up(lineage, window_size_of_interest=200):
    """ Gets rid of the remaining discrepancies """
    
    def recursive_rise_exclusion(lin, totals, rises):
        diff = np.log(lin['length'].values[:-1]) - np.log(
            lin['length'].values[1:])  # What is the difference between the natural logs of two consecutive length measurements across one unit of time?
        
        new_rises = np.where(diff <= -.5)[0]  # Here is where we check if it went up dramatically
        
        new_new = []
        for rise in new_rises:  # for each rise, check if it is followed by an outlier in the next point or not, in which case delete the point before the rise
            # window of interest to see if it is an outlier in that window
            if len(lin) <= window_size_of_interest:  # The raw lineage is too small
                woi = [0, len(lin)]
            elif rise <= window_size_of_interest:  # The raw lineage is big enough but rise is too close to the start of the lineage
                woi = [0, rise + (window_size_of_interest - rise)]
            elif len(lin) - rise <= window_size_of_interest:  # The raw lineage is big enough but rise is too close to the end of the lineage
                woi = [len(lin) - window_size_of_interest, len(lin)]
            else:  # No complications with window size
                woi = [rise - window_size_of_interest, rise + window_size_of_interest]
            
            if lin['length'].iloc[rise + 1] > (lin['length'].iloc[woi[0]:woi[1]].mean() + 2.5 * lin['length'].std()):  # If the next point is abnormally large
                new_new.append(rise + 1)  # replace that point with the one after it, ie. the outlier
            else:
                new_new.append(rise)
        
        new_rises = np.array(new_new)
        
        # if len(new_new) > 0:
        #     plt.plot(lin['length'].values)
        #     plt.scatter(new_rises, lin['length'].values[new_rises], color='brown')
        #     # plt.xlim([2500, 3200])
        #     plt.tight_layout()
        #     plt.show()
        #     plt.close()
        
        rises = np.append(rises, [int(l + np.sum([l >= old for old in totals])) for l in new_rises]).flatten()  # Re-align these point to the raw_data
        
        if len(new_rises) > 0:  # Say we find some abnormal rise
            # Fix the lineage
            new_lin = lin[~lin.index.isin(new_rises)].reset_index(drop=True)
            
            # add to the total indices that are excluded
            new_totals = np.append(totals,
                                   [int(l + np.sum([l >= old for old in totals])) for l in new_rises]).flatten()  # Add the rises indices to the total indices we will ignore for analysis
            
            # Recursion
            lin, totals, rises = recursive_rise_exclusion(new_lin, new_totals, rises)
        
        return [lin, totals, rises]
    
    """ Get rid of the non-positive lengths first """
    total_nans = np.array([])  # Where we will keep all the indices that will be taken out of the corrected lineage
    non_positives = np.array(lineage[lineage['length'] <= 0].index)  # Where the length is non-positive, which is impossible
    
    if len(non_positives) > 0:  # Notify the user and take the non-positives out
        # print('NOTE: This lineage has a length lower than or equal to 0.')
        lineage = lineage[~lineage.index.isin(non_positives)].reset_index(drop=True)  # Goes without saying but there was an instance of this in the Wang Data
    
    total_nans = np.append(total_nans, non_positives).flatten()  # Add the non-positives to the total indices we will ignore for analysis
    
    lineage = lineage.reset_index(drop=True)  # for good-practice
    
    """ Get rid of the single-point singularities first """
    diff = np.log(lineage['length'].values[:-1]) - np.log(
        lineage['length'].values[1:])  # What is the difference between the natural logs of two consecutive length measurements across one unit of time?
    
    # Take out all the singularities
    straight_up = np.where(diff <= -.4)[0]  # These are the points whose next point in time shoot up abnormally
    straight_down = np.where(diff >= .4)[0]  # These are the points whose next point in time fall down abnormally
    
    singularities = np.array([int(down) for up in straight_up for down in straight_down if (down - up == 1)])
    # singularities = np.array([int(down) for up in straight_up for down in straight_down for sep in np.arange(1, 6) if (down - up == sep)])  # Singularities that fall
    singularities = np.append(singularities, np.array([int(up) for up in straight_up for down in straight_down if (up - down == 1)])).flatten()  # Singularities that rise
    # singularities = np.append(singularities, np.array([int(up) + 1 for up in straight_up for down in straight_down if (up - down == 2)])).flatten()  # Singularities that rise
    
    if len(singularities) > 0:  # Notify the user and take the non-positives out
        # print('NOTE: This lineage has singularities that either rise or fall abnormally rapidly.')
        lineage = lineage[~lineage.index.isin(singularities)].reset_index(drop=True)  # Goes without saying but there was an instance of this in the Wang Data
    
    singularities = np.array([int(l + np.sum([l >= old for old in total_nans])) for l in singularities])
    
    total_nans = np.append(total_nans, singularities).flatten()  # Add the singularities to the total indices we will ignore for analysis
    
    """ Get rid of the remaining singularities recursively """
    new_lineage, new_total_nans, failures = recursive_rise_exclusion(lineage, total_nans, rises=np.array([]))
    
    assert len(lineage) >= len(new_lineage)
    assert len(total_nans) <= len(new_total_nans)
    assert (len(singularities) + len(non_positives) + len(failures)) == len(new_total_nans)
    
    return [new_lineage, new_total_nans, failures, singularities, non_positives]


"""  """


def check_the_division(args, lineages=[], raw_lineages=[], raw_indices=[], pu=[]):
    # lineage, cycle_variables_lineage, start_indices, end_indices, raw_lineage, non_positives, singularities, rises, raw_start, raw_end, step_size, total_nans
    
    if args['data_origin'] == '8-31-16 Continue':
        step_size = 6 / 60
    elif args['data_origin'] == 'MG1655_inLB_LongTraces':
        step_size = 5 / 60
    elif args['data_origin'] == 'Maryam_LongTraces':
        step_size = 3 / 60
    elif args['data_origin'] in tanouchi_datasets:
        step_size = 1 / 60
    elif args['data_origin'] in wang_datasets:
        step_size = 1 / 60
    elif args['data_origin'] in sm_datasets:
        step_size = 3 / 60
    else:
        raise IOError('This code is not meant to run the data inputted. Please label the data and put it in as an if-statement.')
    
    if len(raw_lineages) == 0:
        raw_lineages = pd.read_csv(os.path.dirname(os.path.dirname(args['processed_data'])) + '/raw_data_all_in_one.csv')
    if len(raw_indices) == 0:
        raw_indices = pd.read_csv(args['processed_data'] + 'raw_indices_processing.csv')
    if len(pu) == 0:
        pu = pd.read_csv(args['processed_data'] + 'z_score_under_3/physical_units_without_outliers.csv')
    if len(lineages) == 0:
        lineages = raw_lineages.lineage_ID.unique()
    
    seaborn_preamble()
    fig, ax = plt.subplots(tight_layout=True, figsize=[13, 6])
    for lin_id in lineages:
        rl = raw_lineages[raw_lineages['lineage_ID'] == lin_id].copy().sort_values('time').reset_index(drop=True)
        ri = raw_indices[raw_indices['lineage_ID'] == lin_id].copy().sort_values('value').reset_index(drop=True)
        cycle_variables_lineage = pu[pu['lineage_ID'] == lin_id].copy().reset_index(drop=True)
        
        non_positives = ri[ri['type'] == 'non_positives']['value'].values.copy()
        singularities = ri[ri['type'] == 'singularities']['value'].values.copy()
        rises = ri[ri['type'] == 'rises']['value'].values.copy()
        start_indices = ri[ri['type'] == 'start']['value'].values.copy()
        end_indices = ri[ri['type'] == 'end']['value'].values.copy()
        
        plt.plot(rl['length'].values)
        if len(non_positives) > 0:
            print(non_positives)
            print(rl['length'].iloc[non_positives].values)
            plt.scatter(non_positives, rl['length'].iloc[non_positives].values, label='non_positives', marker='o')
        if len(singularities) > 0:
            print(singularities)
            print(rl['length'].iloc[singularities].values)
            plt.scatter(singularities, rl['length'].iloc[singularities].values, label='singularities', marker='v')
        if len(rises) > 0:
            print(rises)
            print(rl['length'].iloc[rises].values)
            plt.scatter(rises, rl['length'].iloc[rises].values, label='failures', marker='^')
        
        for start, end, count in zip(start_indices, end_indices, cycle_variables_lineage.generation.unique()):
            gt = cycle_variables_lineage[cycle_variables_lineage['generation'] == count]['generationtime'].values[0]
            lb = cycle_variables_lineage[cycle_variables_lineage['generation'] == count]['length_birth'].values[0]
            gr = cycle_variables_lineage[cycle_variables_lineage['generation'] == count]['growth_rate'].values[0]
            
            if np.isnan([gt, lb, gr]).any():
                continue
            
            # print([type(l) for l in [gt, gr, lb]])
            #
            # print(step_size)
            #
            # print(type(start), start)
            
            line = lb * np.exp(np.linspace(0, gt, int(np.round(gt * (1 / step_size), 0))) * gr)
            
            # print([type(l) for l in line])
            
            # print(start_indices)
            # print(start, start + len(line))
            # print('8'*9)
            # print(line)
            
            plt.plot(np.arange(start, start + len(line)), line, color='orange', ls='--')  # , label='cycle variable fit'
        
        plt.scatter(start_indices, rl['length'].iloc[start_indices].values, color='green', label='start indices')
        plt.scatter(end_indices, rl['length'].iloc[end_indices].values, color='red', label='end indices')
        plt.legend()
        # plt.tight_layout()
        plt.show()
        plt.close()


""" Create the csv files for physical, trace-centered, and trap-centered units for MM data of a certain type """


def main_mm(args):
    print('type of MM data we are processing:', args['data_origin'])
    
    # Get the file paths to import the data from. 
    # We have to specify the data type, which is different for each experiment.
    file_paths = glob.glob(args['raw_data'] + '/*')
    print(file_paths)
    
    if args['data_origin'] == '20090529_E_coli_Br_SJ119_Wang2010':  # Take out the trajectories that cannot be used
        
        wang_br_sj119_200090529_reject = np.array([
            29, 33, 38, 50, 68, 76, 83, 96, 132, 138, 141, 145, 155, 158, 181, 198, 208, 220, 228, 233, 237, 240, 254, 268, 270, 276, 277, 281, 296, 299
        ]) - 1  # 172, 21, 104, 213
        
        # wang_br_sj119_200090529_cut = {104: 190}
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_br_sj119_200090529_reject)]]
    
    if args['data_origin'] == '20090930_E_coli_MG1655_lexA3_Wang2010':
        wang_MG1655_lexA3_20090930_reject = np.array([
            34, 40, 42, 116
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_lexA3_20090930_reject)]]
    
    if args['data_origin'] == '20090923_E_coli_MG1655_lexA3_Wang2010':
        wang_MG1655_lexA3_20090923_reject = np.array([
            40, 133, 138
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_lexA3_20090923_reject)]]
    
    if args['data_origin'] == '20090922_E_coli_MG1655_lexA3_Wang2010':
        wang_MG1655_lexA3_20090922_reject = np.array([
            4, 40, 133, 138
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_lexA3_20090922_reject)]]
    
    if args['data_origin'] == '20090210_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090210_reject = np.array([
            6, 10
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090210_reject)]]
    
    if args['data_origin'] == '20090129_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090129_reject = np.array([
            6, 10
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090129_reject)]]
    
    if args['data_origin'] == '20090702_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090702_reject = np.array([
            43, 75, 86
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090702_reject)]]
    
    if args['data_origin'] == '20090131_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090131_reject = np.array([
            1, 20, 23, 41, 47, 52, 57
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090131_reject)]]
    
    if args['data_origin'] == '20090525_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090131_reject = np.array([
            89, 105, 112, 127, 233, 250, 318
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090131_reject)]]
    
    if args['data_origin'] == '20090512_E_coli_MG1655_(CGSC_6300)_Wang2010':
        wang_MG1655_CGSC_6300_20090512_reject = np.array([
            54, 64, 135, 146
        ]) - 1
        
        file_paths = np.array(file_paths)[[r for r in np.arange(len(file_paths), dtype=int) if (r not in wang_MG1655_CGSC_6300_20090512_reject)]]
    
    # Create the dataframe for our variables
    cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    with_outliers_cycle_variables = pd.DataFrame(columns=phenotypic_variables + ['lineage_ID', 'generation'])
    
    # The dataframe for our raw data lineages and their indices
    raw_data = pd.DataFrame(columns=['time', 'length', 'lineage_ID', 'filename'])
    raw_indices = pd.DataFrame(columns=['value', 'type', 'lineage_ID'])
    
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
    
    correct_timestamp = ['pos1-1-daughter', 'pos17-1', 'pos17-3', 'pos20', 'pos17-2', 'pos0-1',
                         'pos16-3', 'pos16-2', 'pos16-1', 'pos18-3', 'pos18-2', 'Pos9-1', 'Pos10',
                         'pos19-3', 'Pos9-2', 'pos19-2', 'pos15']
    
    offset = 0  # In case we can't use some files we want the lineage IDs to be in integer order
    
    jump = 0
    
    # load first sheet of each Excel-File, fill internal data structure
    for count, file in enumerate(file_paths[jump:]):
        
        filename = file.split('/')[-1].split('.')[0]
        extension = file.split('/')[-1].split('.')[1]
        
        # Tells us the trap ID and the source (filename)
        print(count + 1 + jump, '/', str(len(file_paths)), ':', filename)
        
        # check=False
        
        # Import the data in a file into a pandas dataframe with the correct extension and pandas function 
        if args['data_origin'] == 'MG1655_inLB_LongTraces':
            step_size = 5 / 60
            # The only file that has an extra column
            if file == 'pos4-4':
                print('In this particular case the lineage divides after the first time-point and it has an extra column.')
                # Extra column when importing
                lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
                # lineage divides after the first time-point
                lineage = lineage.iloc[1:]
            else:  # All the rest don't have these problems
                lineage = pd.read_csv(file, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            lineage['time'] = lineage['time'] * (5 / 3)  # Incorrect labelling
        elif args['data_origin'] == 'Maryam_LongTraces':
            step_size = 3 / 60
            # This is because some the data is in .xls format while others are in .csv
            if extension == 'csv':
                lineage = pd.read_csv(file, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
            elif extension == 'xls':
                lineage = pd.read_excel(file, names=['time', 'length'])[['time', 'length']].dropna(axis=0)
            else:
                raise IOError('For MaryamLongTraces dataset non-xls/csv files have not been inspected!')
        elif args['data_origin'] in tanouchi_datasets:
            lineage = pd.read_csv(file, delimiter=',', names=['time', 'division_flag', 'length', 'fluor', 'avg_fluor'])
            lineage['time'] = (lineage['time'] - 1) / 60  # Because we map the index to the correct time-step-size which is 1 minute
            step_size = 1 / 60  # one-minute measurements!
        # elif args['data_origin'] == '8-31-16 Continue':
        #     lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        #
        #     # print(lineage['time'])
        #     # print(lineage['time'].iloc[0])
        #     # print('_'*100)
        #
        #     if .05 == lineage['time'].iloc[1]:
        #         # Wrong labeling of time
        #         lineage['time'] = lineage['time'] * 2
        #     step_size = 6 / 60
        elif args['data_origin'] == 'lambda_LB':
            # There are quite a lot of files with an extra column at the beginning
            if filename in extra_column:
                lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            elif filename == 'pos15':
                print('This one is special.')
                lineage = pd.read_csv(file, delimiter='\t', names=['_', 'time', 'length', 'something similar to length', 'something protein', 'other protein', '__', '___', '___1'])[
                    ['time', 'length']]
            else:
                lineage = pd.read_csv(file, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
            
            if filename in correct_timestamp:
                lineage['time'] = lineage['time'] * 2
            
            step_size = max(np.unique(np.diff(lineage['time'])), key=list(np.diff(lineage['time'])).count)
            
            print('step size', step_size)
        elif args['data_origin'] in wang_datasets:
            lineage = pd.read_csv(file, delimiter=' ')  # names=['time', 'division_flag', 'length', 'width', 'area', 'yfp_intensity', 'CMx', 'CMy']
            
            # Sometimes the time index is called time and sometimes its called index
            if 'time' in lineage.columns:
                lineage = lineage.rename(columns={'division': 'division_flag'})[['time', 'division_flag', 'length']]
            elif 'index' in lineage.columns:
                lineage = lineage.rename(columns={'index': 'time', 'division': 'division_flag'})[['time', 'division_flag', 'length']]
            
            # We have to make sure that the indices are monotonically increasing by 1 in order to trust the time axis
            if len(np.unique(np.diff(lineage['time']))) != 1 or np.unique(np.diff(lineage['time']))[0] != 1:
                # print(np.unique(np.diff(lineage['time'])))
                # raise IOError('time given in Wang data is not monotonous and increasing.')
                print('the time given in the data file is not increasing always by 1')  # so we do not know how to measure time for this lineage we will not use it.')
                check = True
                # continue
            
            lineage['time'] = (lineage['time']) / 60  # Because we map the index to the correct time-step-size which is 1 minute
            
            lineage['length'] = (lineage['length']) * 0.0645  # Convert it from pixel length to micrometers
            
            step_size = 1 / 60  # one-minute measurements!
            
            # plt.plot(lineage['length'].values)
            # plt.show()
            # plt.close()
            
            if (args['data_origin'] == '20090702_E_coli_MG1655_(CGSC_6300)_Wang2010'):  # Manual cleanup so that the clean_up algorithm works properly
                if (filename == 'xy08_ch2_cell0'):
                    lineage = lineage.drop([1175, 1176, 1177, 1178], axis=0)
                elif (filename == 'xy04_ch0_cell0'):
                    lineage = lineage.drop([995], axis=0)
                elif (filename == 'xy10_ch6_cell0'):
                    lineage = lineage.drop([1034, 1035], axis=0)
                elif (filename == 'xy09_ch7_cell0'):
                    lineage = lineage.iloc[:2628]
                elif (filename == 'xy03_ch16_cell0'):
                    lineage = lineage.drop([701, 702], axis=0)
                    lineage = lineage.iloc[:2165]
                elif (filename == 'xy05_ch10_cell0'):
                    lineage = lineage.drop([2326, 2327, 2328], axis=0)
                elif (filename == 'xy10_ch10_cell0'):
                    lineage = lineage.iloc[:2451]
                    lineage = lineage.drop([1603, 1604], axis=0)
                    lineage = lineage.drop([1588, 1589, 1590], axis=0)
                    lineage = lineage.drop([1002, 1003], axis=0)
                    lineage = lineage.drop(np.arange(2015, 2111), axis=0)
                elif (filename == 'xy06_ch3_cell0'):
                    lineage = lineage.drop([1467, 1468, 1547, 1548], axis=0)
                elif (filename == 'xy09_ch10_cell0'):
                    lineage = lineage.drop([1427, 1428], axis=0)
                elif (filename == 'xy01_ch4_cell0'):
                    lineage = lineage.drop(np.arange(1565, 1571), axis=0)
                elif (filename == 'xy01_ch2_cell0'):
                    lineage = lineage.drop(np.arange(1039, 1045), axis=0)
                    lineage = lineage.drop([544, 545, 546], axis=0)
                elif (filename == 'xy02_ch9_cell0'):
                    lineage = lineage.drop([336, 337, 338], axis=0)
                elif (filename == 'xy04_ch6_cell0'):
                    lineage = lineage.iloc[:2628]
                elif (filename == 'xy09_ch17_cell0'):
                    lineage = lineage.iloc[:2367]
                elif (filename == 'xy06_ch2_cell0'):
                    lineage = lineage.drop([772, 773], axis=0)
                elif (filename == 'xy06_ch8_cell0'):
                    lineage = lineage.iloc[:2743]
                elif (filename == 'xy09_ch16_cell0'):
                    lineage = lineage.iloc[:2073]
                elif (filename == 'xy06_ch4_cell0'):
                    lineage = lineage.drop([2312, 2310], axis=0)
                elif (filename == 'xy10_ch2_cell0'):
                    lineage = lineage.iloc[:1834]
                elif (filename == 'xy09_ch8_cell0'):
                    lineage = lineage.iloc[:2042]
                elif (filename == 'xy09_ch3_cell0'):
                    lineage = lineage.drop([1565, 1566], axis=0)
                    lineage = lineage.drop(np.arange(1481, 1486), axis=0)
                    lineage = lineage.drop([1420, 1421, 1442, 1443], axis=0)
                    lineage = lineage.drop([350, 351, 352, 353], axis=0)
                elif (filename == 'xy10_ch4_cell0'):
                    lineage = lineage.iloc[:2686]
                elif (filename == 'xy09_ch13_cell0'):
                    lineage = lineage.drop([2324, 2325], axis=0)
                elif (filename == 'xy10_ch15_cell0'):
                    lineage = lineage.drop(np.arange(1987, 1995), axis=0)
                elif (filename == 'xy04_ch15_cell0'):
                    lineage = lineage.drop([962], axis=0)
                elif (filename == 'xy10_ch9_cell0'):
                    lineage = lineage.drop([2159, 2160], axis=0)
                if (filename == 'xy05_ch12_cell0'):
                    lineage = lineage.drop(np.arange(611, 616), axis=0)
                if (filename == 'xy09_ch14_cell0'):
                    lineage = lineage.drop(np.arange(1307, 1311), axis=0)
                if (filename == 'xy04_ch13_cell0'):
                    lineage = lineage.drop(np.arange(2604, 2625), axis=0)
                    lineage = lineage.drop(np.arange(1487, 1538), axis=0)
                if (filename == 'xy06_ch6_cell0'):
                    lineage = lineage.drop(np.arange(1121, 1123), axis=0)
                if (filename == 'xy03_ch3_cell0'):
                    lineage = lineage.drop([1160, 1161], axis=0)
                if (filename == 'xy08_ch14_cell0'):
                    lineage = lineage.drop(np.arange(710, 719), axis=0)
            if (args['data_origin'] == '20090131_E_coli_MG1655_(CGSC_6300)_Wang2010'):
                if (filename == 'xy13c1_ch0_cell0'):
                    lineage = lineage.iloc[:2992]
                elif (filename == 'xy07c1_ch0_cell0'):
                    lineage = lineage.drop(np.arange(2794, 2944), axis=0)
            if (args['data_origin'] == '20090525_E_coli_MG1655_(CGSC_6300)_Wang2010'):
                if (filename == 'xy13_ch7_cell0_YFP0002'):
                    lineage = lineage.iloc[:319]
                elif (filename == 'xy09_ch5_cell0_YFP0002'):
                    lineage = lineage.iloc[:592]
                elif (filename == 'xy11_ch0_cell0_YFP0001'):
                    lineage = lineage.drop(np.arange(1094, 1294), axis=0)
                elif (filename == 'xy07_ch0_cell0_YFP0002'):
                    lineage = lineage.iloc[:969]
                elif (filename == 'xy02_ch14_cell0_YFP0001'):
                    lineage = lineage.iloc[:2118]
                elif (filename == 'xy08_ch5_cell0_YFP0001'):
                    lineage = lineage.iloc[:2993]
                elif (filename == 'xy02_ch3_cell0_YFP0002'):
                    lineage = lineage.drop(np.arange(884, 891), axis=0)
                    lineage = lineage.drop(np.arange(783, 800), axis=0)
                elif (filename == 'xy06_ch4_cell0_YFP0002'):
                    lineage = lineage.drop(np.arange(73, 78), axis=0)
                elif (filename == 'xy14_ch3_cell0_YFP0002'):
                    lineage = lineage.iloc[:277]
                elif (filename == 'xy06_ch5_cell0_YFP0001'):
                    lineage = lineage.drop(np.arange(295, 297), axis=0)
            if (args['data_origin'] == '20090512_E_coli_MG1655_(CGSC_6300)_Wang2010'):
                if (filename == 'xy05_ch8_cell0'):
                    lineage = lineage.drop(np.arange(1122, 1126), axis=0)
                if (filename == 'xy11_ch9_cell0'):
                    lineage = lineage.drop([540, 541], axis=0)
                if (filename == 'xy04_ch11_cell0'):
                    lineage = lineage.drop([1377, 1378], axis=0)
                if (filename == 'xy05_ch19_cell0'):
                    lineage = lineage.drop([699, 700], axis=0)
                if (filename == 'xy09_ch8_cell0'):
                    lineage = lineage.drop(np.arange(1515, 1520), axis=0)
                if (filename == 'xy05_ch12_cell0'):
                    lineage = lineage.drop(np.arange(1322, 1325), axis=0)
        
        else:
            raise IOError('This code is not meant to run the data inputted. Please label the data and put it in as an if-statement.')
        
        # Check if time is going backwards and discard if true
        if not all(x < y for x, y in zip(lineage['time'].values[:-1], lineage['time'].values[1:])):
            print(filename, ': Time is going backwards. We cannot use this data.')
            
            # reset the lineage_ID
            offset += 1
            continue
        
        # elif args['data_origin'] == 'LAC_M9':
        #     # Simple :)
        #     raw_lineage = pd.read_csv(filename, delimiter='\t', names=['time', 'length', 'something similar to length', 'something protein', 'other protein'])[['time', 'length']]
        # Make sure we have the measurement time step-size in hours and that it is the same across all rows
        # If not, do not use the trace (Too much of a headache for the phenotypic variable linear regression).
        # step_sizes = (raw_lineage['time'].iloc[1:].values - raw_lineage['time'].iloc[:-1].values).round(2)
        # if not np.all(step_sizes == step_sizes[0]):
        #     print(filename, ' has steps that are not the same size ', np.unique(step_sizes), ' and we are not using this data then')  # , are we are not gonna use these...')
        #     continue
        #     # exit()
        # if args['data_origin'] == 'lambda_LB' and .05 == step_sizes[0]:
        #     pass
        #     step_size = .05
        # elif args['data_origin'] == 'lambda_LB' and .1 == step_sizes[0]:
        #     print('step size == .1', filename, np.unique(step_sizes))
        #     raw_lineage['time'] = raw_lineage['time'] / 2
        #     step_size = .05
        # elif args['data_origin'] == 'lambda_LB':
        #     print('step size == .06 repeating', filename, np.unique(step_sizes))
        #     raw_lineage['time'] = (raw_lineage['time'] * 15) / 20
        #     step_size = .05
        
        lineage['filename'] = filename  # for reference
        lineage['lineage_ID'] = count + 1 - offset  # Add the lineage ID
        lineage['step_size'] = step_size
        
        raw_lineage = lineage.copy()  # Copy the lineage before we edit it
        
        # add it to the total data
        raw_data = raw_data.append(lineage, ignore_index=True)
        
        lineage, total_nans, rises, singularities, non_positives = clean_up(lineage)
        # total_nans, rises, singularities, non_positives = [], [], [], []
        
        # These datasets have a division flag
        # if args['data_origin'] in tanouchi_datasets:  # + wang_datasets
        #     start_indices = np.array(raw_lineage[raw_lineage['division_flag'] == 1].index)
        #     end_indices = np.append(start_indices[1:] - 1, len(raw_lineage) - 1)
        #     # exit()
        # else:
        #     # Figure out the indices for the division events
        #     start_indices, end_indices = get_division_indices(raw_lineage['length'].values)
        
        try:
            start_indices, end_indices = get_division_indices(lineage['length'].values)  # Figure out the indices for the division events
            # print(start_indices)
            # print(total_nans)
            # start_indices = start_indices[~np.isin(start_indices, total_nans)]
            # end_indices = start_indices[~np.isin(start_indices, total_nans)]
        except:
            continue
        
        # plt.plot(np.arange(len(raw_lineage['length']), dtype=int), raw_lineage['length'], color='blue', label='raw_lineage')
        # # plt.plot(np.arange(start_indices[0], len(interpolation)+start_indices[0]), interpolation, color='orange', label='cycle variable fit', ls='--')
        # plt.scatter(start_indices, raw_lineage['length'].iloc[start_indices], color='green', label='start indices')
        # plt.scatter(end_indices, raw_lineage['length'].iloc[end_indices], color='red', label='end indices')
        # plt.title(filename)
        # plt.xlabel('index')
        # plt.ylabel(r'length ($\mu$m)')
        # plt.yscale('log')
        # plt.tight_layout()
        # plt.show()
        # plt.close()
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage, with_outliers, start_indices, \
        end_indices, outsider_start_indices, outsider_end_indices = linear_regression(lineage, raw_lineage, start_indices, end_indices, int(count + 1 - offset), fit_the_lengths=True)
        
        outsider_start_indices = [int(start + np.sum([start >= l for l in total_nans])) for start in outsider_start_indices]
        outsider_end_indices = [int(end + np.sum([end >= l for l in total_nans])) for end in outsider_end_indices]
        
        raw_start = [int(
            start + np.sum([start >= l for l in rises]) + np.sum([start >= l for l in singularities]) + np.sum([start >= l for l in non_positives])
        ) for start in start_indices]
        raw_end = [int(
            end + np.sum([end >= l for l in rises]) + np.sum([end >= l for l in singularities]) + np.sum([end >= l for l in non_positives])
        ) for end in end_indices]
        
        # raw_start = [int(start + np.sum([start >= l for l in total_nans])) for start in start_indices]
        # raw_end = [int(end + np.sum([end >= l for l in total_nans])) for end in end_indices]
        
        to_append = {}
        blah = 0
        for value, label in zip([outsider_start_indices, outsider_end_indices, raw_start, raw_end, rises, singularities, non_positives],
                                ['outsider', 'outsider', 'start', 'end', 'rises', 'singularities', 'non_positives']):
            for v in value:
                to_append.update({  # append the start
                    blah: {
                        'value': v,
                        'type': label,
                        'lineage_ID': count + 1 - offset
                    }
                })
                
                blah += 1
        
        raw_indices = raw_indices.append(pd.DataFrame.from_dict(to_append, "index"), ignore_index=True)
        
        # if np.array([0 != l for l in [len(total_nans), len(non_positives), len(singularities), len(rises)]]).any():
        #     print('{} number of ignored points\n{} number of non-positives\n{} number of singularities\n{} number of long-failures'.format(len(total_nans), len(non_positives), len(singularities),
        #                                                                                                                                    len(rises)))
        
        # if check:
        #     check_the_division(args, lineages=[count + 1 - offset], raw_lineages=raw_lineage, raw_indices=pd.DataFrame.from_dict(to_append, "index"), pu=cycle_variables_lineage)
        # check_the_division(args, lineages=[count + 1 - offset], raw_lineages=raw_lineage, raw_indices=pd.DataFrame.from_dict(to_append, "index"), pu=cycle_variables_lineage)
        
        # append the cycle variables to the
        cycle_variables = cycle_variables.append(cycle_variables_lineage, ignore_index=True)
        with_outliers_cycle_variables = with_outliers_cycle_variables.append(with_outliers, ignore_index=True)
        
        print('-' * 200)
    
    print('processed data:\n', cycle_variables)
    print('cleaned raw data:\n', raw_data)
    
    for variable in phenotypic_variables:
        if variable == 'division_ratio':
            print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - len(cycle_variables['lineage_ID'].unique())))
        else:
            print(variable, cycle_variables[variable].count() / len(cycle_variables[variable]))
    
    # Save the raw data
    raw_data.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(os.path.dirname(os.path.dirname(args['raw_data'])) + '/raw_data_all_in_one.csv', index=False)
    
    # Save the ones with the z score under 3
    without_outliers = args['processed_data'] + 'z_score_under_3'
    create_folder(without_outliers)
    cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(without_outliers + '/physical_units_without_outliers.csv', index=False)
    minusing(cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        without_outliers + '/trace_centered_without_outliers.csv',
        index=False)
    
    # Save the ones with outliers still
    with_outliers_cycle_variables.reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(args['processed_data'] + '/physical_units.csv', index=False)
    minusing(with_outliers_cycle_variables.reset_index(drop=True), phenotypic_variables).reset_index(drop=True).sort_values(['lineage_ID', 'generation']).to_csv(
        args['processed_data'] + '/trace_centered.csv',
        index=False)
    
    # Save the raw indices that are supposed to be highlighted
    raw_indices.sort_values(['lineage_ID']).reset_index(drop=True).to_csv(args['processed_data'] + '/raw_indices_processing.csv', index=False)


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
    raw_indices = pd.DataFrame(columns=['value', 'type', 'lineage_ID'])
    
    # The files that have the RawData
    files = glob.glob(args['raw_data'] + '/*.xls')
    
    # load first sheet of each Excel-File, fill rawdata dataframe
    for count, file in enumerate(files):
        # print(count, file.split('/')[-1], sep=': ')
        
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
    
    step_size = .05
    
    for lineage_id in raw_data.lineage_ID.unique():
        print('Lineage ID:', lineage_id)
        
        # Get the lineage
        lineage = raw_data[raw_data['lineage_ID'] == lineage_id].copy()
        
        raw_lineage = lineage.copy()  # Copy the lineage before we edit it
        
        lineage, total_nans, rises, singularities, non_positives = clean_up(lineage)
        
        # Figure out the indices for the division events
        start_indices, end_indices = get_division_indices(lineage['length'].values)
        
        # add the cycle variables to the overall dataframe
        cycle_variables_lineage, with_outliers, start_indices, \
        end_indices, outsider_start_indices, outsider_end_indices = linear_regression(lineage, raw_lineage, start_indices, end_indices, int(lineage_id),
                                                                                      fit_the_lengths=True)
        
        # plt.plot(lineage['length'].values)
        # plt.scatter(start_indices, lineage['length'].iloc[start_indices])
        # plt.scatter(end_indices, lineage['length'].iloc[end_indices])
        # plt.show()
        # plt.close()
        
        # Add the SM categorical variables
        cycle_variables_lineage['trap_ID'] = lineage['trap_ID'].unique()[0]
        cycle_variables_lineage['trace'] = lineage['trace'].unique()[0]
        cycle_variables_lineage['dataset'] = lineage['dataset'].unique()[0]
        
        with_outliers['trap_ID'] = lineage['trap_ID'].unique()[0]
        with_outliers['trace'] = lineage['trace'].unique()[0]
        with_outliers['dataset'] = lineage['dataset'].unique()[0]
        
        outsider_start_indices = [int(start + np.sum([start >= l for l in total_nans])) for start in outsider_start_indices]
        outsider_end_indices = [int(end + np.sum([end >= l for l in total_nans])) for end in outsider_end_indices]
        
        raw_start = [int(start + np.sum([start >= l for l in total_nans])) for start in start_indices]
        raw_end = [int(end + np.sum([end >= l for l in total_nans])) for end in end_indices]
        
        to_append = {}
        blah = 0
        for value, label in zip([outsider_start_indices, outsider_end_indices, raw_start, raw_end, rises, singularities, non_positives],
                                ['outsider', 'outsider', 'start', 'end', 'rises', 'singularities', 'non_positives']):
            for v in value:
                to_append.update({  # append the start
                    blah: {
                        'value': v,
                        'type': label,
                        'lineage_ID': lineage_id
                    }
                })
                
                blah += 1
        
        raw_indices = raw_indices.append(pd.DataFrame.from_dict(to_append, "index"), ignore_index=True)
        
        if np.array([l != 0 for l in [len(total_nans), len(non_positives), len(singularities), len(rises)]]).any():
            print('{} number of ignored points\n{} number of non-positives\n{} number of singularities\n{} number of long-failures'.format(len(total_nans), len(non_positives), len(singularities),
                                                                                                                                           len(rises)))
        
        # Append the cycle variables to the processed dataframe
        cycle_variables = cycle_variables.append(cycle_variables_lineage[order], ignore_index=True)
        with_outliers_cycle_variables = with_outliers_cycle_variables.append(with_outliers[order], ignore_index=True)
    
    print('processed data:\n', cycle_variables)
    
    print('cleaned raw data:\n', raw_data)
    
    # Check how much NaNs were introduced because of the zscore < 3 condition on one of the dataframes (no outliers)
    for variable in phenotypic_variables:
        if variable in ['division_ratio', 'div_and_fold']:  # , 'fold_then_div'
            # This is because, by definition, the first two variables have 1 NaN value at the first generation, while the third variable has 1 NaN at the end of each lineage
            print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - (1 * len(cycle_variables['lineage_ID'].unique()))))
        # elif variable in ['div_then_fold']:
        #     # This is because, by definition, div then fold variable has two NaNs at the first two generations of each lineage
        #     print(variable, cycle_variables[variable].count() / (len(cycle_variables[variable]) - (2 * len(cycle_variables['lineage_ID'].unique()))))
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
    
    # Save the raw indices that are supposed to be highlighted
    raw_indices.reset_index(drop=True).sort_values(['lineage_ID']).to_csv(args['processed_data'] + '/raw_indices_processing.csv', index=False)


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
    for data_origin in ['Maryam_LongTraces']:  # input_args.dataset_names[21:]:
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
        
        # check_the_division(args, lineages=[275, 273, 266, 263])


if __name__ == '__main__':
    main()
