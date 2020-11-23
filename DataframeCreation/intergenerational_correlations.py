#!/usr/bin/env bash

import numpy as np
import pandas as pd
from os import mkdir
import os
import argparse
import sys
from scipy.stats import pearsonr
import time

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets, boot_pearson


""" Creates the mom and daughter young and old synced csv files for the correlation df creation """


def intergenerational_dataframe_creations_SM(info_columns, info, how_many_separations, directory):
    old_kind_specific_df = pd.DataFrame(columns=info_columns)
    young_kind_specific_df = pd.DataFrame(columns=info_columns)
    
    # 1 separation is mother daughter, 2 is grandmother granddaughter, etc...
    for separation in range(0, how_many_separations + 1):
        # the dataframe for the older and the younger generations
        older_df = pd.DataFrame(columns=info_columns)
        younger_df = pd.DataFrame(columns=info_columns)
        
        # loop over every key in both S and NS datasets
        for dataset in np.unique(info['dataset']):
            for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
                both_traj_df = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id)].sort_values(['trace', 'generation'])
                a_trace = both_traj_df[both_traj_df['trace'] == 'A'].reset_index(drop=True)
                b_trace = both_traj_df[both_traj_df['trace'] == 'B'].reset_index(drop=True)
                
                # decides which generation to put so we don't have any conflicting or non-true inter-generational pairs
                older_mask = [ind for ind in range(len(a_trace) - separation)] + [len(a_trace) + ind for ind in range(len(b_trace) - separation)]
                younger_mask = [separation + ind for ind in range(len(a_trace) - separation)] + [separation + len(a_trace) + ind for ind in
                                                                                                 range(len(b_trace) - separation)]
                
                # add this trap's mother and daughter cells to the df
                older_df = pd.concat([older_df, both_traj_df.iloc[older_mask]], axis=0, join='inner').reset_index(drop=True)
                younger_df = pd.concat([younger_df, both_traj_df.iloc[younger_mask]], axis=0, join='inner').reset_index(drop=True)
        
        young_kind_specific_df = young_kind_specific_df.append(younger_df, ignore_index=True)
        old_kind_specific_df = old_kind_specific_df.append(older_df, ignore_index=True)
        
        # # save them to csb
        # older_df.to_csv('{}/Older Dataframe of distance {}.csv'.format(directory, separation), index=False)
        # younger_df.to_csv('{}/Younger Dataframe of distance {}.csv'.format(directory, separation), index=False)
        #
        # print('distance {}'.format(separation))

    young_kind_specific_df.to_csv('{}/young_inter_correlations.csv'.format(args.save_folder), index=False)
    old_kind_specific_df.to_csv('{}/old_inter_correlations.csv'.format(args.save_folder), index=False)


""" Creates the mom and daughter young and old synced csv files for the correlation df creation """


def intergenerational_creations(info_columns, info, how_many_separations, label):
    def get_df(separation, young_param, old_param, dataset, corr_df, n_boots, older_df, younger_df):
        to_add = dict(zip(list(corr_df.columns), [[] for _ in range(len(corr_df.columns))]))
            
        # # read the csvs
        # older_df = pd.read_csv('{}/Older Dataframe of distance {}.csv'.format(pu_dir, separation))
        # younger_df = pd.read_csv('{}/Younger Dataframe of distance {}.csv'.format(pu_dir, separation))
        
        if n_boots == 0:
            to_add['distance'] = np.append(to_add['distance'], separation)
            to_add['correlations'] = np.append(to_add['correlations'], pearsonr(older_df[old_param], younger_df[young_param])[0])
            to_add['dataset'] = np.append(to_add['dataset'], dataset)
            to_add['old_var'] = np.append(to_add['old_var'], old_param)
            to_add['young_var'] = np.append(to_add['young_var'], young_param)
            
            # convert the arrays into numpy arrays
            to_add['distance'] = np.array(to_add['distance']).flatten()
            to_add['correlations'] = np.array(to_add['correlations']).flatten()
            to_add['dataset'] = np.array(to_add['dataset']).flatten()
            to_add['old_var'] = np.array(to_add['old_var']).flatten()
            to_add['young_var'] = np.array(to_add['young_var']).flatten()
        else:
            
            # Put them into a dictionary
            to_add['distance'] = np.append(to_add['distance'], [separation for _ in range(n_boots)])
            to_add['correlations'] = np.append(to_add['correlations'], boot_pearson(older_df[old_param], younger_df[young_param], n_boots))
            to_add['dataset'] = np.append(to_add['dataset'], [dataset for _ in range(n_boots)])
            to_add['old_var'] = np.append(to_add['old_var'], [old_param for _ in range(n_boots)])
            to_add['young_var'] = np.append(to_add['young_var'], [young_param for _ in range(n_boots)])
            
            # convert the arrays into numpy arrays
            to_add['distance'] = np.array(to_add['distance']).flatten()
            to_add['correlations'] = np.array(to_add['correlations']).flatten()
            to_add['dataset'] = np.array(to_add['dataset']).flatten()
            to_add['old_var'] = np.array(to_add['old_var']).flatten()
            to_add['young_var'] = np.array(to_add['young_var']).flatten()
        
        # print(to_add)
        # exit()
        
        # return the new dataframe added to the old one
        return pd.concat([corr_df, pd.DataFrame(to_add)], axis=0)
    
    old_kind_specific_df = pd.DataFrame(columns=info_columns)
    young_kind_specific_df = pd.DataFrame(columns=info_columns)
    
    output_df = pd.DataFrame(columns=['distance', 'correlations', 'dataset', 'old_var', 'young_var'])
    
    # 1 separation is mother daughter, 2 is grandmother granddaughter, etc...
    for separation in range(1, how_many_separations + 1):
        # the dataframe for the older and the younger generations
        older_df = pd.DataFrame(columns=info_columns)
        younger_df = pd.DataFrame(columns=info_columns)
        
        # loop over every key in both S and NS datasets
        for dataset in np.unique(info['dataset']):
            print(dataset)
            for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
                both_traj_df = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id)]
                a_trace = both_traj_df[both_traj_df['trace'] == 'A'].sort_values(['trace', 'generation']).reset_index(drop=True)
                b_trace = both_traj_df[both_traj_df['trace'] == 'B'].sort_values(['trace', 'generation']).reset_index(drop=True)
                
                # decides which generation to put so we don't have any conflicting or non-true inter-generational pairs
                older_mask = [ind for ind in range(len(a_trace) - separation)] + [len(a_trace) + ind for ind in range(len(b_trace) - separation)]
                younger_mask = [separation + ind for ind in range(len(a_trace) - separation)] + [separation + len(a_trace) + ind for ind in
                                                                                                 range(len(b_trace) - separation)]

                # add this trap's mother and daughter cells to the df
                older_df = pd.concat([older_df, both_traj_df.iloc[older_mask]], axis=0, join='inner').reset_index(drop=True)
                younger_df = pd.concat([younger_df, both_traj_df.iloc[younger_mask]], axis=0, join='inner').reset_index(drop=True)
        
        for old_param in phenotypic_variables:
            # print(old_param)
            for young_param, ind in zip(phenotypic_variables, range(len(phenotypic_variables))):
                # print(older_df[old_param], younger_df[young_param])
                # print(older_df, younger_df)
                # exit()
                output_df = get_df(separation, young_param, old_param, label, output_df, n_boots, older_df, younger_df)
                # total_df = total_df.append(output_df, ignore_index=True)
                # output_df = get_df(separation, young_param, old_param, 'trace centered', output_df, n_boots, older_df, younger_df)
                # print(output_df)
                # exit()

        # output_df = output_df.append()
        
    return output_df


""" Creates a csv file which contains bootstrapped correlations of paramter pairings """


def correlation_dataframe_creation_SM():
    def get_df(how_many_separations, young_param, old_param, dataset, output_df, n_boots, older_df, younger_df):
        to_add = dict(zip(list(output_df.columns), [[] for _ in range(len(output_df.columns))]))
        
        # 1 separation is mother daughter, 2 is grandmother granddaughter, etc...
        for separation in range(0, how_many_separations + 1):
            
            # # read the csvs
            # older_df = pd.read_csv('{}/Older Dataframe of distance {}.csv'.format(pu_dir, separation))
            # younger_df = pd.read_csv('{}/Younger Dataframe of distance {}.csv'.format(pu_dir, separation))
            
            if n_boots == 0:
                to_add['distance'] = np.append(to_add['distance'], separation)
                to_add['correlations'] = np.append(to_add['correlations'], pearsonr(older_df[old_param], younger_df[young_param])[0])
                to_add['dataset'] = np.append(to_add['dataset'], dataset)
                to_add['old_var'] = np.append(to_add['old_var'], old_param)
                to_add['young_var'] = np.append(to_add['young_var'], young_param)

                # convert the arrays into numpy arrays
                to_add['distance'] = np.array(to_add['distance']).flatten()
                to_add['correlations'] = np.array(to_add['correlations']).flatten()
                to_add['dataset'] = np.array(to_add['dataset']).flatten()
                to_add['old_var'] = np.array(to_add['old_var']).flatten()
                to_add['young_var'] = np.array(to_add['young_var']).flatten()
            else:
            
                # Put them into a dictionary
                to_add['distance'] = np.append(to_add['distance'], [separation for _ in range(n_boots)])
                to_add['correlations'] = np.append(to_add['correlations'], boot_pearson(older_df[old_param], younger_df[young_param], n_boots))
                to_add['dataset'] = np.append(to_add['dataset'], [dataset for _ in range(n_boots)])
                to_add['old_var'] = np.append(to_add['old_var'], [old_param for _ in range(n_boots)])
                to_add['young_var'] = np.append(to_add['young_var'], [young_param for _ in range(n_boots)])
                
                # convert the arrays into numpy arrays
                to_add['distance'] = np.array(to_add['distance']).flatten()
                to_add['correlations'] = np.array(to_add['correlations']).flatten()
                to_add['dataset'] = np.array(to_add['dataset']).flatten()
                to_add['old_var'] = np.array(to_add['old_var']).flatten()
                to_add['young_var'] = np.array(to_add['young_var']).flatten()
            
            # print(to_add)
            # exit()
        
        # return the new dataframe added to the old one
        return pd.concat([output_df, pd.DataFrame(to_add)], axis=0)
    
    young_df = pd.read_csv('{}/young_inter_correlations.csv'.format(args.save_folder))
    old_df = pd.read_csv('{}/old_inter_correlations.csv'.format(args.save_folder))
    
    output_df = pd.DataFrame(columns=['distance', 'correlations', 'dataset', 'old_var', 'young_var'])
    
    # Put 1 because we know that both their mothers are the same
    for old_param in phenotypic_variables:
        print(old_param)
        for young_param in phenotypic_variables:
            # output_df = pd.DataFrame(columns=['distance', 'correlations', 'dataset'])
            
            output_df = get_df(how_many_separations, young_param, old_param, 'physical units', output_df, n_boots, old_df, young_df)
            # total_df = total_df.append(output_df, ignore_index=True)
            output_df = get_df(how_many_separations, young_param, old_param, 'trace centered', output_df, n_boots, old_df, young_df)
            # total_df = total_df.append(output_df, ignore_index=True)
            
            # output_df.to_csv('{}/{}, {}.csv'.format(corr_dir, old_param, young_param), index=False)
            # print('Saved {}, {}.csv'.format(old_param, young_param))
    output_df.to_csv('{}/total_inter_correlations.csv'.format(args.save_folder), index=False)


""" script that creates the correlation dataframes """


start_time = time.time()

    
parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered.csv')
parser.add_argument('-bs', '--bootstraps', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)

args = parser.parse_args()

create_folder(args.save_folder)

# import the labeled measured bacteria in trace-centered units
physical_units = pd.read_csv(args.pu)
trace_centered = pd.read_csv(args.tc)

# Here we create the directory where we will save the sorted dataframes
pu_dir = '{}/old and young dataframes physical units'.format(args.save_folder)
create_folder(pu_dir)

print('physical units')

# # Create the synchronized csv files
# intergenerational_dataframe_creations_SM(physical_units.columns, physical_units, how_many_separations=15, directory=pu_dir)
#
# print("--- %s seconds ---" % (time.time() - start_time))

# Here we create the directory where we will save the sorted dataframes
tc_dir = '{}/old and young dataframes trace centered'.format(args.save_folder)
create_folder(tc_dir)

print('trace-centered')


n_boots = 0

print('pu')
total_inter_df = intergenerational_creations(physical_units.columns, physical_units, how_many_separations=15, label='physical units').reset_index(drop=True)
print('tc')
total_inter_df = total_inter_df.append(intergenerational_creations(trace_centered.columns, trace_centered, how_many_separations=15, label='trace centered'), ignore_index=True).reset_index(drop=True)

total_inter_df.to_csv('{}/total_inter_correlations.csv'.format(args.save_folder), index=False)

# # Create the synchronized csv files
# intergenerational_dataframe_creations_SM(trace_centered.columns, trace_centered, how_many_separations=15, directory=tc_dir)

print("--- %s seconds ---" % (time.time() - start_time))
exit()


# Where we will store our bootstrapped pearson correlations
corr_dir = '{}/correlation dataframes'.format(args.save_folder)
create_folder(corr_dir)

how_many_separations = 15
n_boots = 0

# Now we create the correlation dataframes
correlation_dataframe_creation_SM()

print("--- %s seconds ---" % (time.time() - start_time))

