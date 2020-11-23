#!/usr/bin/env bash

import sys
sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, shuffle_info
import pandas as pd
import numpy as np


""" Creates a dataframe that contains all the kl divergences """


def kl_divergence(info, phenotypic_variables, kind):
    # kl div of two univariate normal distributions
    def kl(m0, m1, c0, c1):
        return np.log(c1 / c0) + (c0 ** 2 + (m0 - m1) ** 2) / (2 * (c1 ** 2)) - .5
    
    # Initialize the dataframe where we will keep the kl-divergences
    kl_df = pd.DataFrame(columns=['value', 'variable', 'trap_ID', 'trace', 'kind'])
    
    for param in phenotypic_variables:
        # The population average
        pop_mean = info[param].mean()
        pop_cov = info[param].std()
        
        for dataset in np.unique(info['dataset']):
            
            for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
                
                for trace in ['A', 'B']:
                    
                    # mean and covariance of the
                    trace_mean = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id) & (info['trace'] == trace)][param].mean()
                    trace_cov = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id) & (info['trace'] == trace)][param].std()
                    
                    # Because KL-Divergence is not symmetric we do two types and see how similar they are...
                    kldiv1 = kl(pop_mean, trace_mean, pop_cov, trace_cov)
                    kldiv2 = kl(trace_mean, pop_mean, trace_cov, pop_cov)
                    
                    # add it to the dataframe where we keep it symmetrized
                    kl_df = kl_df.append(
                        {'value': kldiv1, 'variable': param, 'trap_ID': trap_id, 'trace': trace, 'kind': kind},
                        ignore_index=True)
                    kl_df = kl_df.append(
                        {'value': kldiv2, 'variable': param, 'trap_ID': trap_id, 'trace': trace, 'kind': kind},
                        ignore_index=True)
    
    return kl_df


""" Returns dataframe with same size containing the time-averages of each phenotypic variable instead of the local value """


def get_time_averages_df(info, phenotypic_variables):
    # We keep the trap means here
    means_df = pd.DataFrame(columns=['dataset', 'trap_ID', 'trace', 'max_gen', 'generation'] + phenotypic_variables)
    
    # specify a lineage
    for dataset in ['SP', 'NC']:
        for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
            for trace in ['A', 'B']:
                
                # the values of the lineage we get from physical units
                lineage = info[(info['trap_ID'] == trap_id) & (info['dataset'] == dataset) & (info['trace'] == trace)].copy()
                
                # add its time-average
                to_add = {
                    'dataset': [dataset for _ in np.arange(len(lineage))], 'trap_ID': [trap_id for _ in np.arange(len(lineage))], 'trace': [trace for _ in np.arange(len(lineage))],
                    'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
                }
                to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables})
                to_add = pd.DataFrame(to_add)
                means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    
    assert len(info) == len(means_df)
    
    return means_df


""" Creates a dataframe with the ergodicity breaking parameter of each phenotypic variable """


def ergodicity_breaking_parameter(df, phenotypic_variables, kind, n_boots=0):
    # Initialize where we will put the bootstrapped ergodicity breaking variable
    eb_df = pd.DataFrame(columns=['variable', 'kind', 'value'])
    
    # get the dataframes where in each entry, instead of the generation specific value, there is the time-average
    time_averages = get_time_averages_df(df, phenotypic_variables)
    
    # to normalize the different phenotypic_variables
    pop_var = df.var()
    
    # bootstrap this ergodicity breaking parameter
    if n_boots != 0:
        for _ in np.arange(n_boots):
            # Bootstrapping the indices has the lineage length inequality taken into account
            indices = time_averages.sample(frac=1, replace=True).index
            
            # Get the variance of the time-averages for both kinds of lineages
            variance = time_averages[phenotypic_variables].loc[indices].var() / pop_var
            
            # add them both to the dataframe where we save it all
            for param in phenotypic_variables:
                eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    else:
        # Get the variance of the time-averages for both kinds of lineages
        variance = time_averages[phenotypic_variables].var() / pop_var
    
        # add them both to the dataframe where we save it all
        for param in phenotypic_variables:
            eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    
    return [time_averages, eb_df]


import argparse


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')


args = parser.parse_args()


create_folder(args.save_folder)


# import the labeled measured bacteria in physical units
info = pd.read_csv('{}/physical_units.csv'.format(args.save_folder))


print('Population physical_units lineages')
# Create lineages sampled from a population distribution
shuffled = shuffle_info(info)
shuffled.to_csv('{}/PopulationLineages.csv'.format(args.save_folder), index=False)


print('ergodicity breaking and time-averages')
# get the bootstrapped EB variable for both kinds of lineages
time_averages_trace, eb_df = ergodicity_breaking_parameter(info, phenotypic_variables, kind='Trace')
_, eb_df_pop = ergodicity_breaking_parameter(shuffled, phenotypic_variables, kind='Population')
eb_df = eb_df.append(eb_df_pop, ignore_index=True).reset_index(drop=True)


# save it to the right folder
time_averages_trace.to_csv('{}/time_averages.csv'.format(args.save_folder), index=False)
eb_df.to_csv('{}/ergodicity_breaking_parameter.csv'.format(args.save_folder), index=False)


print('kl_divergences')
# Put in the kl divergences for each variable for each type of lineage
kl_df = kl_divergence(info, phenotypic_variables, 'Trace')
kl_df = kl_df.append(kl_divergence(shuffled, phenotypic_variables, 'Population'), ignore_index=True).reset_index(drop=True)


# save the kl_df dataframe
kl_df.to_csv('{}/kullback_leibler_divergences.csv'.format(args.save_folder), index=False)
