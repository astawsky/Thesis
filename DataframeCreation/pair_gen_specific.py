#!/usr/bin/env bash

import sys
import os
import pandas as pd
import numpy as np
import argparse
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets

""" Create the correlation between all pair bacteria across lineages in a dataset """
def generation_specific(physical_units, trace_means_df, number_of_generations, bootstraps):
    
    """ Get the covariance and correlation and error from theoretical equality """
    
    # dataframe that we will save
    df = pd.DataFrame(
        columns=['variable_a', 'variable_b', 'AB_parameter_pair', 'covariance', 'correlation', 'kind',
                 'cov difference between data and model' 'dataset', 'generation'])
    
    # dataframe that we will save
    bs = pd.DataFrame(
        columns=['variable_a', 'variable_b', 'AB_parameter_pair', 'covariance', 'correlation', 'kind',
                 'cov difference between data and model', 'dataset', 'generation'])
    
    for dataset in datasets:
        # go over all symmetric pairings
        
        print(dataset)
        
        repeats = []
        for variable_a in phenotypic_variables:
            for variable_b in phenotypic_variables:
                if variable_b not in repeats:
                # if variable_b == variable_a:
                    
                    # An approximation so that we normalize all generations by the same std and more robust that taking it at each generation
                    a_std = np.mean([np.std(physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == gen)][variable_a]) for gen in np.arange(number_of_generations)])
                    b_std = np.mean([np.std(physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == gen)][variable_b]) for gen in np.arange(number_of_generations)])
                    
                    for generation in np.arange(number_of_generations):
                        # The dataset mask
                        pu_a = physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[True, True])[
                            variable_a].reset_index(drop=True)
                        pu_b = physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[False, True])[
                            variable_b].reset_index(drop=True)
                        
                        ta_a = trace_means_df[(trace_means_df['dataset'] == dataset) & (trace_means_df['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[True, True])[
                            variable_a].reset_index(drop=True)
                        ta_b = trace_means_df[(trace_means_df['dataset'] == dataset) & (trace_means_df['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[False, True])[
                            variable_b].reset_index(drop=True)
                        
                        # do it without bootstrap so we have an empirical measurement of the time-averages
                        pu_cov = np.cov(pu_a, pu_b)[0, 1]
                        ta_cov = np.cov(ta_a, ta_b)[0, 1]
                        tc_cov = pu_cov - ta_cov
                        actual_cov = ta_cov + tc_cov
                        
                        # Calculate the "correlations"
                        pu_corr = pu_cov / (a_std * b_std)
                        ta_corr = ta_cov / (a_std * b_std)
                        tc_corr = tc_cov / (a_std * b_std)
                        actual_corr = actual_cov / (a_std * b_std)
                        
                        # Add them to the dataframe we will export
                        for kind, cov, corr in [
                            ('physical_units', pu_cov, pu_corr),
                            ('trace_centered', tc_cov, tc_corr),
                            ('time_averages', ta_cov, ta_corr),
                            ('model', actual_cov, actual_corr)
                        ]:
                            
                            # the labels
                            if variable_a == variable_b:
                                # It looks nicer for the paper
                                label = '{}'.format(symbols[kind][variable_a])
                            else:
                                label = '{}, {}'.format(symbols[kind][variable_a], symbols[kind][variable_b])
                            
                            # save it to the dataframe
                            df = df.append({
                                'variable_a': variable_a,
                                'variable_b': variable_b,
                                'AB_parameter_pair': label,
                                'covariance': cov,
                                'correlation': corr,
                                'kind': kind,
                                'cov difference between data and model': actual_cov - pu_cov,
                                'corr difference between data and model': actual_corr - pu_corr,
                                'dataset': dataset,
                                'generation': generation + 1
                            }, ignore_index=True)
                        
                        if not bootstraps:
                            pass
                        else:
                            # do the bootstrap
                            for _ in np.arange(bootstraps):
                                indices = \
                                    physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[True, True])[
                                        variable_a].reset_index(
                                        drop=True).sample(frac=1, replace=True).index
                                
                                # The dataset mask
                                pu_a = physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[True, True])[
                                    variable_a].reset_index(drop=True).iloc[
                                    indices]
                                pu_b = \
                                    physical_units[(physical_units['dataset'] == dataset) & (physical_units['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[False, True])[
                                        variable_b].reset_index(drop=True).iloc[
                                        indices]
                                
                                ta_a = trace_means_df[(trace_means_df['dataset'] == dataset) & (trace_means_df['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[True, True])[
                                    variable_a].reset_index(drop=True).iloc[
                                    indices]
                                ta_b = \
                                    trace_means_df[(trace_means_df['dataset'] == dataset) & (trace_means_df['generation'] == generation)].sort_values(['trace', 'trap_ID'], ascending=[False, True])[
                                        variable_b].reset_index(drop=True).iloc[
                                        indices]
                                
                                # do it without bootstrap so we have an empirical measurement of the time-averages
                                pu_cov = np.cov(pu_a, pu_b)[0, 1]  # np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(pu_a, pu_b)])
                                ta_cov = np.cov(ta_a, ta_b)[0, 1]  # np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(ta_a, ta_b)])
                                tc_cov = pu_cov - ta_cov  # np.mean([a * b for a, b in zip(tc_a, tc_b)])
                                actual_cov = pu_cov - ta_cov
                                
                                pu_corr = pu_cov / (a_std * b_std)
                                tc_corr = tc_cov / (a_std * b_std)
                                ta_corr = ta_cov / (a_std * b_std)
                                actual_corr = actual_cov / (a_std * b_std)
                                
                                for kind, cov, cor in [
                                    ('physical_units', pu_cov, pu_corr),
                                    ('trace_centered', tc_cov, tc_corr),
                                    ('time_averages', ta_cov, ta_corr),
                                    ('model', actual_cov, actual_corr)
                                ]:
                                    # the labels
                                    if variable_a == variable_b:
                                        # It looks nicer for the paper
                                        label = '{}'.format(symbols[kind][variable_a])
                                    else:
                                        label = '{}, {}'.format(symbols[kind][variable_a], symbols[kind][variable_b])
                                    
                                    # save it to the dataframe
                                    bs = bs.append({
                                        'variable_a': variable_a,
                                        'variable_b': variable_b,
                                        'AB_parameter_pair': label,
                                        'covariance': cov,
                                        'correlation': corr,
                                        'kind': kind,
                                        'cov difference between data and model': actual_cov - pu_cov,
                                        'corr difference between data and model': actual_corr - pu_corr,
                                        'dataset': dataset,
                                        'generation': generation + 1
                                    }, ignore_index=True)
            
            # append the parameter so we do not have any repetitions
            repeats.append(variable_a)
    
    return [df, bs]


def main(args):
    # get the three datasets we need to compute the covariance and correlations
    physical_units, trace_means_df = pd.read_csv('{}/{}'.format(args.save_folder, args.puc)), pd.read_csv('{}/{}'.format(args.save_folder, args.tac))
    
    pair_correlation, bootstrapped = generation_specific(physical_units, trace_means_df, args.gen_specific_amount, args.bs)
    
    # save it to the Data folder
    pair_correlation.to_csv('{}/over_lineages.csv'.format(args.save_folder), index=False)
    if not args.bs:
        pass
    else:
        bootstrapped.to_csv('{}/over_lineages_{}_bootstraps.csv'.format(args.save_folder, args.bs), index=False)


# parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
# parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
# parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units_with_control.csv')
# parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages_with_control.csv')
# parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered_with_control.csv')
# parser.add_argument('-bs', '--bs', metavar='', type=int, help='How many bootstraps per covariance should be done?',
#                     required=False, default=0)
#
# args = parser.parse_args()
#
# create_folder(args.save_folder)
#
# # get the three datasets we need to compute the covariance and correlations
# physical_units, trace_centered, trace_means_df = pd.read_csv(args.pu), pd.read_csv(args.tc), pd.read_csv(args.ta)
#
# number_of_generations = 6
#
# pair_correlation, bootstrapped = generation_specific(bootstraps=args.bs)
#
# # save it to the Data folder
# pair_correlation.to_csv('{}/over_lineages.csv'.format(args.save_folder), index=False)
# if not args.bs:
#     pass
# else:
#     bootstrapped.to_csv('{}/over_lineages_{}_bootstraps.csv'.format(args.save_folder, args.bs), index=False)
