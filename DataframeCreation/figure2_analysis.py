#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, shuffle_info, shuffle_lineage_generations, limit_lineage_length, trace_center_a_dataframe
import pandas as pd
import numpy as np

""" adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """


def expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df, label):
    import matplotlib.pyplot as plt
    
    # create the info dataframe with values being the cumulative mean of the lineage
    for dataset in np.unique(df['dataset']):
        
        for trap_id in np.unique(df[df['dataset'] == dataset]['trap_ID']):
            # Go through both traces in an SM trap
            for trace_id in ['A', 'B']:
                # specify the trace
                trace = df[(df['trap_ID'] == int(trap_id)) & (df['trace'] == trace_id) & (df['dataset'] == dataset)].copy()
                
                # add its time-average up until and including this generation
                to_add = pd.DataFrame.from_dict({
                    'label': [label for _ in range(len(trace))],
                    'trap_ID': [int(trap_id) for _ in range(len(trace))],
                    'trace': [trace_id for _ in range(len(trace))],
                    'generation': [generation + 1 for generation in trace['generation']]
                }).reset_index(drop=True)

                # plt.axhline(0, color='black')
                # plt.plot(trace.sort_values('generation')['generationtime'].expanding().mean().values, label='expanding mean')
                # plt.plot(trace.sort_values('generation')['generationtime'].cumsum().values, label='cumsum')
                # plt.legend()
                # plt.show()
                # plt.close()
                # exit()
                
                to_add_expanding_mean = pd.concat([trace.sort_values('generation')[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
                to_add_cum_sum = pd.concat([trace.sort_values('generation')[phenotypic_variables].cumsum().reset_index(drop=True), to_add], axis=1)
                expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
                cum_sum_df = cum_sum_df.append(to_add_cum_sum, ignore_index=True).reset_index(drop=True)
                
                if to_add_expanding_mean.isnull().values.any():
                    print('expanding mean')
                    print(dataset, trap_id, trace)
                    print(to_add_expanding_mean)
                    input()
                if cum_sum_df.isnull().values.any():
                    print('expanding mean')
                    print(dataset, trap_id, trace)
                    print(cum_sum_df)
                    input()
                
    assert not expanding_mean.isnull().values.any()
    assert not cum_sum_df.isnull().values.any()
    
    # Calculate the cv and var over all lineages in the dataset
    for param in phenotypic_variables:
        for generation in np.arange(1, df['generation'].max()+1):
            time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)][param]
            different_walks = cum_sum_df[(cum_sum_df['label'] == label) & (cum_sum_df['generation'] == generation)][param]
            
            if time_averages.isnull().values.any():
                print('expanding mean')
                print(time_averages)
                input()
            if different_walks.isnull().values.any():
                print('expanding mean')
                print(different_walks)
                input()
            
            # print(time_averages)
            # print(different_walks)
            # print(time_averages.isnull().values.any())
            # print(different_walks.isnull().values.any())
            # exit()
            
            # adding it to the dataframe of all coefficients of variations
            variance_of_expanding_mean = variance_of_expanding_mean.append(
                {
                    'label': label,
                    'generation': generation,
                    'cv': time_averages.std() / time_averages.mean(),
                    'var': time_averages.var(),
                    'param': param
                }, ignore_index=True
            )
            variance_of_cum_sum_df = variance_of_cum_sum_df.append(
                {
                    'label': label,
                    'generation': generation,
                    'var': different_walks.var(),
                    'param': param
                }, ignore_index=True
            )
    
    return [expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df]


def main(args):
    # # import/create the trace lineages
    # physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    #
    # physical_units = limit_lineage_length(physical_units, min_gens=50)
    # trace_centered = trace_center_a_dataframe(physical_units)
    
    # print(physical_units.isnull().any())
    # print(trace_centered.isnull().any())
    # print(len(physical_units), len(trace_centered))
    # exit()
    
    # import/create the trace lineages
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))

    # import/create the trace-centered lineages
    trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tc))
    
    # import/create the population lineages
    try:
        population_sampled = pd.read_csv('{}/{}'.format(args.save_folder, args.population_sampled))
        exit()
    except:
        print('creating population sampled')
        population_sampled = shuffle_info(physical_units, MM=args.MM)
        population_sampled.to_csv('{}/{}'.format(args.save_folder, args.population_sampled), index=False)
    
    # import/create the shuffled generations lineages
    try:
        shuffled_generations = pd.read_csv('{}/{}'.format(args.save_folder, args.shuffled))
        exit()
    except:
        print('creating shuffled generations')
        shuffled_generations = shuffle_lineage_generations(physical_units)
        shuffled_generations.to_csv('{}/{}'.format(args.save_folder, args.shuffled), index=False)
    
    # The generation-shuffled trace-centered dataframe
    try:
        shuffled_tc = pd.read_csv('{}/{}'.format(args.save_folder, args.tc_shuffled))
        exit()
    except:
        print('creating shuffled tc')
        shuffled_tc = shuffle_lineage_generations(trace_centered)
        shuffled_tc.to_csv('{}/{}'.format(args.save_folder, args.tc_shuffled), index=False)
    
    # We keep the trap means here
    expanding_mean = pd.DataFrame(columns=['label', 'trap_ID', 'trace', 'generation'] + phenotypic_variables)
    
    # Keep the cv per lineage length here here
    variance_of_expanding_mean = pd.DataFrame(columns=['label', 'generation', 'cv', 'var', 'param'])
    
    # We keep the trap means here
    cum_sum_df = pd.DataFrame(columns=['label', 'trap_ID', 'trace', 'generation'] + phenotypic_variables)
    
    # Keep the cv per lineage length here here
    variance_of_cum_sum_df = pd.DataFrame(columns=['label', 'generation', 'var', 'param'])
    
    # Calculate the cv and TA per lineage length
    for kind, df in zip(['Shuffled', 'Trace', 'Population', 'Trace-Centered', 'Shuffled TC'], [shuffled_generations, physical_units, population_sampled, trace_centered, shuffled_tc]):
        # kind = 'Trace-Centered'
        # df = trace_centered
        print(kind)
        expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df = expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, variance_of_expanding_mean,
                                                                                                                             cum_sum_df,
                                                                                                                             variance_of_cum_sum_df, kind)

    # Save the csv file
    expanding_mean.to_csv('{}/{}'.format(args.save_folder, args.cta), index=False)

    # Save the csv file
    variance_of_expanding_mean.to_csv('{}/{}'.format(args.save_folder, args.vcta), index=False)

    # Save the csv file
    cum_sum_df.to_csv('{}/{}'.format(args.save_folder, args.cum_sum), index=False)

    # Save the csv file
    variance_of_cum_sum_df.to_csv('{}/{}'.format(args.save_folder, args.v_cum_sum), index=False)
