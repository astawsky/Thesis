#!/usr/bin/env bash

import sys
import os
import pandas as pd
import numpy as np
import argparse

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets

""" creates a dataframe that contains all the kl divergences """


def kl_divergence():
    def kl(m0, m1, c0, c1):
        return .5 * (np.trace(np.dot(np.linalg.inv(c1), c0)) + np.dot((m1 - m0).T, np.dot(np.linalg.inv(c1), (m1 - m0))) - len(m0) + np.log((np.linalg.det(c1) / np.linalg.det(c0))))
    
    # where we keep the kldivergences
    kl_df = pd.DataFrame(columns=['value', 'variable', 'trap_ID', 'dataset'])
    
    for dataset in np.unique(info['dataset']):
        for variable in phenotypic_variables:
            for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
                
                # mean and covariance of the pairs as Gaussian Parameter approximations
                a_mean = info[(info['trap_ID'] == trap_id) & (info['trace'] == 'A') & (info['dataset'] == dataset)][[variable]].mean()
                a_cov = info[(info['trap_ID'] == trap_id) & (info['trace'] == 'A') & (info['dataset'] == dataset)][[variable]].cov()
                b_mean = info[(info['trap_ID'] == trap_id) & (info['trace'] == 'B') & (info['dataset'] == dataset)][[variable]].mean()
                b_cov = info[(info['trap_ID'] == trap_id) & (info['trace'] == 'B') & (info['dataset'] == dataset)][[variable]].cov()
                
                # Because KL-Divergence is not symmetric we do two types and see how similar they are...
                # We compare the A and B trace, not the population-ensemble
                kldiv1 = kl(a_mean, b_mean, a_cov, b_cov)
                kldiv2 = kl(b_mean, a_mean, b_cov, a_cov)
                
                # add it to the dataframe where we keep it all
                kl_df = kl_df.append(
                    {'value': kldiv1, 'variable': variable, 'trap_ID': trap_id, 'dataset': dataset},
                    ignore_index=True)
                kl_df = kl_df.append(
                    {'value': kldiv2, 'variable': variable, 'trap_ID': trap_id, 'dataset': dataset},
                    ignore_index=True)
    
    return kl_df


""" Create the correlation between all pair bacteria across lineages in a dataset """


def pair_bacteria(bootstraps=False):
    # get the three datasets we need to compute the covariance and correlations
    physical_units, trace_centered, trace_means_df = pd.read_csv(args.pu), pd.read_csv(args.tc), pd.read_csv(args.ta)
    
    # dataframe that we will save
    df = pd.DataFrame(
        columns=['variable_a', 'variable_b', 'AB_parameter_pair', 'covariance', 'correlation', 'kind', 'correlation error of average from theory', 'dataset'])
    
    # dataframe that we will save
    bs = pd.DataFrame(
        columns=['variable_a', 'variable_b', 'AB_parameter_pair', 'covariance', 'correlation', 'kind', 'correlation error of average from theory', 'dataset'])
    
    for dataset in datasets:
        # go over all symmetric pairings
        
        print(dataset)
        
        repeats = []
        for variable_a in phenotypic_variables:
            print(variable_a)
            for variable_b in phenotypic_variables:
                # if variable_b not in repeats:
                if variable_b == variable_a:
                    
                    # The dataset mask
                    pu_a = physical_units[(physical_units['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(drop=True)
                    pu_b = physical_units[(physical_units['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(drop=True)
                    
                    tc_a = trace_centered[(trace_centered['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(drop=True)
                    tc_b = trace_centered[(trace_centered['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(drop=True)
                    
                    ta_a = trace_means_df[(trace_means_df['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(drop=True)
                    ta_b = trace_means_df[(trace_means_df['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(drop=True)
                    
                    a_mean = np.mean(ta_a)
                    b_mean = np.mean(ta_b)
                    
                    a_std = np.std(pu_a)
                    b_std = np.std(pu_b)
                    
                    # important from the theory
                    assert np.mean(ta_a).round(8) == np.mean(ta_b).round(8) == np.mean(pu_a).round(8) == np.mean(pu_b).round(8)
                    
                    # do it without bootstrap so we have an empirical measurement of the time-averages
                    pu_cov = np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(pu_a, pu_b)])  # np.cov(pu_a, pu_b)[0, 1]  #
                    ta_cov = np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(ta_a, ta_b)])  # np.cov(ta_a, ta_b)[0, 1]  #
                    tc_cov = np.mean([a * b for a, b in zip(tc_a, tc_b)])  # np.cov(tc_a, tc_b)[0, 1]  #
                    actual_cov = ta_cov + tc_cov
                    
                    pu_corr = pu_cov / (a_std * b_std)
                    tc_corr = tc_cov / (a_std * b_std)
                    ta_corr = ta_cov / (a_std * b_std)
                    actual_corr = actual_cov / (a_std * b_std)
                    
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
                            'dataset': dataset
                        }, ignore_index=True)
                    
                    if not bootstraps:
                        pass
                    else:
                        # do the bootstrap
                        for _ in np.arange(bootstraps):
                            indices = physical_units[(physical_units['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(
                                drop=True).sample(
                                frac=1, replace=True).index
                            
                            # The dataset mask
                            pu_a = \
                                physical_units[(physical_units['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(
                                    drop=True).loc[
                                    indices]
                            pu_b = \
                                physical_units[(physical_units['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(
                                    drop=True).loc[
                                    indices]
                            
                            tc_a = \
                                trace_centered[(trace_centered['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(
                                    drop=True).loc[
                                    indices]
                            tc_b = \
                                trace_centered[(trace_centered['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(
                                    drop=True).loc[
                                    indices]
                            
                            ta_a = \
                                trace_means_df[(trace_means_df['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, True, True])[variable_a].reset_index(
                                    drop=True).loc[
                                    indices]
                            ta_b = \
                                trace_means_df[(trace_means_df['dataset'] == dataset)].sort_values(['generation', 'trace', 'trap_ID'], ascending=[True, False, True])[variable_b].reset_index(
                                    drop=True).loc[
                                    indices]
                            
                            # do it without bootstrap so we have an empirical measurement of the time-averages
                            pu_cov = np.cov(pu_a, pu_b)[0, 1]  # np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(pu_a, pu_b)])
                            ta_cov = np.cov(ta_a, ta_b)[0, 1]  # np.mean([(a - a_mean) * (b - b_mean) for a, b in zip(ta_a, ta_b)])
                            tc_cov = np.cov(tc_a, tc_b)[0, 1]  # np.mean([a * b for a, b in zip(tc_a, tc_b)])
                            actual_cov = ta_cov + tc_cov
                            
                            pu_corr = pu_cov / (a_std * b_std)
                            tc_corr = tc_cov / (a_std * b_std)
                            ta_corr = ta_cov / (a_std * b_std)
                            actual_corr = actual_cov / (a_std * b_std)
                            
                            for kind, cov, corr, percentage in [
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
                                    'dataset': dataset
                                }, ignore_index=True)
            
            # append the parameter so we do not have any repetitions
            repeats.append(variable_a)
    
    return [df, bs]


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units_with_control.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages_with_control.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered_with_control.csv')
parser.add_argument('-bs', '--bootstraps', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)

args = parser.parse_args()

create_folder(args.save_folder)

# import the labeled measured bacteria in physical units
info = pd.read_csv('{}/physical_units_with_control.csv'.format(args.save_folder))

# Put in the kl divergences for each parameter for each type of lineage
kld_df = kl_divergence()

# save the kl_df dataframe
kld_df.to_csv('{}/kullback_leibler_divergences_for_pair_lineages.csv'.format(args.save_folder), index=False)

pair_correlation, bootstrapped = pair_bacteria(bootstraps=args.bootstraps)

# save it to the Data folder
pair_correlation.to_csv('{}/over_lineages_and_time.csv'.format(args.save_folder), index=False)
if not args.bootstraps:
    pass
else:
    bootstrapped.to_csv('{}/over_lineages_and_time_{}_bootstraps.csv'.format(args.save_folder, args.bootstraps), index=False)
