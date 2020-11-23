#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, shuffle_info, shuffle_lineage_generations
import pandas as pd
import numpy as np


""" adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """


def cumulative_mean_and_cv_dataframes(df, phenotypic_variables, trap_means_df, cv_df, label):
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
                
                to_add = pd.concat([trace.sort_values('generation')[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)  # index=np.arange(len(trace), dtype=int)
                trap_means_df = trap_means_df.append(to_add, ignore_index=True).reset_index(drop=True)
                
                # for generation in trace['generation']:
                #     # add its time-average up until and including this generation
                #     to_add = {'label': label, 'trap_ID': trap_id, 'trace': trace_id, 'generation': generation + 1}
                #     to_add.update(dict(zip(phenotypic_variables, trace[(trace['generation'] <= generation)].mean())))
                #     trap_means_df = trap_means_df.append(to_add, ignore_index=True)
    
    # Calculate the cv and var over all lineages in the dataset
    for param in phenotypic_variables:
        print(param)
        for generation in np.arange(1, 51):
            time_averages = trap_means_df[(trap_means_df['label'] == label) & (trap_means_df['generation'] == generation)][param]
            
            # adding it to the dataframe of all coefficients of variations
            cv_df = cv_df.append(
                {
                    'label': label,
                    'generation': generation,
                    'cv': time_averages.std() / time_averages.mean(),
                    'var': time_averages.var(),
                    'param': param
                }, ignore_index=True
            )
    
    return [trap_means_df, cv_df]


import argparse

parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')

args = parser.parse_args()

create_folder(args.save_folder)


pd.set_option("display.max_columns", None)


# import/create the trace lineages
physical_units = pd.read_csv('{}/physical_units.csv'.format(args.save_folder))

# import/create the trace-centered lineages
trace_centered = pd.read_csv('{}/trace_centered.csv'.format(args.save_folder))

# import/create the population lineages
try:
    population_sampled = pd.read_csv('{}/PopulationLineages.csv'.format(args.save_folder))
except:
    population_sampled = shuffle_info(physical_units)
    population_sampled.to_csv('{}/PopulationLineages.csv'.format(args.save_folder), index=False)

# import/create the shuffled generations lineages
try:
    shuffled_generations = pd.read_csv('{}/shuffled_generations.csv'.format(args.save_folder))
except:
    shuffled_generations = shuffle_lineage_generations(physical_units)
    shuffled_generations.to_csv('{}/shuffled_generations.csv'.format(args.save_folder), index=False)


# We keep the trap means here
trap_means_df = pd.DataFrame(columns=['label', 'trap_ID', 'trace', 'generation'] + phenotypic_variables)

# Keep the cv per lineage length here here
cv_df = pd.DataFrame(columns=['label', 'generation', 'cv', 'var', 'param'])

# Calculate the cv and TA per lineage length
print('Generation shuffled')
trap_means_df, cv_df = cumulative_mean_and_cv_dataframes(shuffled_generations, phenotypic_variables, trap_means_df, cv_df, 'Shuffled')
print('-' * 50)
print('Trace')
trap_means_df, cv_df = cumulative_mean_and_cv_dataframes(physical_units, phenotypic_variables, trap_means_df, cv_df, 'Trace')
print('-' * 50)
print('Population')
trap_means_df, cv_df = cumulative_mean_and_cv_dataframes(population_sampled, phenotypic_variables, trap_means_df, cv_df, 'Population')
print('-' * 50)
print('Trace-Centered')
trap_means_df, cv_df = cumulative_mean_and_cv_dataframes(trace_centered, phenotypic_variables, trap_means_df, cv_df, 'Trace-Centered')

# Save the csv file
trap_means_df.to_csv('{}/cumulative_time_averages.csv'.format(args.save_folder), index=False)

# Save the csv file
cv_df.to_csv('{}/variation_of_cumulative_time_averages.csv'.format(args.save_folder), index=False)
