#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables_new, shuffle_info, dataset_names
import pandas as pd
import numpy as np


""" Creates a dataframe that contains all the kl divergences """


def kl_divergence(info, phenotypic_variables_new, kind):
    # kl div of two univariate normal distributions
    def kl(m0, m1, c0, c1):
        return np.log(c1 / c0) + (c0 ** 2 + (m0 - m1) ** 2) / (2 * (c1 ** 2)) - .5

    # Initialize the dataframe where we will keep the kl-divergences
    kl_df = pd.DataFrame(columns=['value', 'variable', 'lineage_ID', 'kind'])

    for param in phenotypic_variables_new:
        # The population average
        pop_mean = info[param].mean()
        pop_cov = info[param].std()

        for lineage_id in info['lineage_ID'].unique():
    
            # mean and covariance of the
            trace_mean = info[(info['lineage_ID'] == lineage_id)][param].mean()
            trace_cov = info[(info['lineage_ID'] == lineage_id)][param].std()
    
            # Because KL-Divergence is not symmetric we do two types and see how similar they are...
            kldiv1 = kl(pop_mean, trace_mean, pop_cov, trace_cov)
            kldiv2 = kl(trace_mean, pop_mean, trace_cov, pop_cov)
    
            # add it to the dataframe where we keep it symmetrized
            kl_df = kl_df.append(
                {'value': kldiv1, 'variable': param, 'lineage_ID': lineage_id, 'kind': kind},
                ignore_index=True)
            kl_df = kl_df.append(
                {'value': kldiv2, 'variable': param, 'lineage_ID': lineage_id, 'kind': kind},
                ignore_index=True)
    
    # if MM:
    #     # Initialize the dataframe where we will keep the kl-divergences
    #     kl_df = pd.DataFrame(columns=['value', 'variable', 'lineage_ID', 'kind'])
    #
    #     for param in phenotypic_variables_new:
    #         # The population average
    #         pop_mean = info[param].mean()
    #         pop_cov = info[param].std()
    #
    #         for trap_id in info['lineage_ID'].unique():
    #
    #             # mean and covariance of the
    #             trace_mean = info[(info['lineage_ID'] == trap_id)][param].mean()
    #             trace_cov = info[(info['lineage_ID'] == trap_id)][param].std()
    #
    #             # Because KL-Divergence is not symmetric we do two types and see how similar they are...
    #             kldiv1 = kl(pop_mean, trace_mean, pop_cov, trace_cov)
    #             kldiv2 = kl(trace_mean, pop_mean, trace_cov, pop_cov)
    #
    #             # add it to the dataframe where we keep it symmetrized
    #             kl_df = kl_df.append(
    #                 {'value': kldiv1, 'variable': param, 'lineage_ID': trap_id, 'kind': kind},
    #                 ignore_index=True)
    #             kl_df = kl_df.append(
    #                 {'value': kldiv2, 'variable': param, 'lineage_ID': trap_id, 'kind': kind},
    #                 ignore_index=True)
    # else:
    #     # Initialize the dataframe where we will keep the kl-divergences
    #     kl_df = pd.DataFrame(columns=['value', 'variable', 'trap_ID', 'trace', 'kind'])
    #
    #     for param in phenotypic_variables_new:
    #         # The population average
    #         pop_mean = info[param].mean()
    #         pop_cov = info[param].std()
    #
    #         for dataset in np.unique(info['dataset']):
    #
    #             for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
    #
    #                 for trace in ['A', 'B']:
    #
    #                     # mean and covariance of the
    #                     trace_mean = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id) & (info['trace'] == trace)][param].mean()
    #                     trace_cov = info[(info['dataset'] == dataset) & (info['trap_ID'] == trap_id) & (info['trace'] == trace)][param].std()
    #
    #                     # Because KL-Divergence is not symmetric we do two types and see how similar they are...
    #                     kldiv1 = kl(pop_mean, trace_mean, pop_cov, trace_cov)
    #                     kldiv2 = kl(trace_mean, pop_mean, trace_cov, pop_cov)
    #
    #                     # add it to the dataframe where we keep it symmetrized
    #                     kl_df = kl_df.append(
    #                         {'value': kldiv1, 'variable': param, 'trap_ID': trap_id, 'trace': trace, 'kind': kind},
    #                         ignore_index=True)
    #                     kl_df = kl_df.append(
    #                         {'value': kldiv2, 'variable': param, 'trap_ID': trap_id, 'trace': trace, 'kind': kind},
    #                         ignore_index=True)
    
    return kl_df


""" Returns dataframe with same size containing the time-averages of each phenotypic variable instead of the local value """


def get_time_averages_df(info, phenotypic_variables_new):  # MM
    # We keep the trap means here
    means_df = pd.DataFrame(columns=['lineage_ID', 'max_gen', 'generation'] + phenotypic_variables_new)
    
    # specify a lineage
    for trap_id in info['lineage_ID'].unique():
        
        # the values of the lineage we get from physical units
        lineage = info[(info['lineage_ID'] == trap_id)].copy()
        
        # add its time-average
        to_add = {
            'lineage_ID': [trap_id for _ in np.arange(len(lineage))],
            'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
        }
        to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables_new})
        to_add = pd.DataFrame(to_add)
        means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    
    assert len(info) == len(means_df)
    
    # if MM:
    #     # We keep the trap means here
    #     means_df = pd.DataFrame(columns=['lineage_ID', 'max_gen', 'generation'] + phenotypic_variables_new)
    #
    #     # specify a lineage
    #     for trap_id in info['lineage_ID'].unique():
    #
    #         # the values of the lineage we get from physical units
    #         lineage = info[(info['lineage_ID'] == trap_id)].copy()
    #
    #         # add its time-average
    #         to_add = {'lineage_ID': [trap_id for _ in np.arange(len(lineage))],
    #             'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
    #         }
    #         to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables_new})
    #         to_add = pd.DataFrame(to_add)
    #         means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    #
    #     assert len(info) == len(means_df)
    # else:
    #     # We keep the trap means here
    #     means_df = pd.DataFrame(columns=['dataset', 'trap_ID', 'trace', 'max_gen', 'generation'] + phenotypic_variables_new)
    #
    #     # specify a lineage
    #     for lineage_id in info.lineage_ID.unique():
    #         # the values of the lineage we get from physical units
    #         lineage = info[(info['lineage_ID'] == lineage_id)].copy()
    #
    #         # dataset = lineage.dataset.unique()[0]
    #         # trap_id = lineage.trap_ID.unique()[0]
    #         # trace = lineage.trace.unique()[0]
    #         #
    #         # # add its time-average
    #         # to_add = {
    #         #     'dataset': [dataset for _ in np.arange(len(lineage))], 'trap_ID': [trap_id for _ in np.arange(len(lineage))], 'trace': [trace for _ in np.arange(len(lineage))],
    #         #     'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
    #         # }
    #
    #         # add its time-average
    #         to_add = {
    #             'dataset': [dataset for _ in np.arange(len(lineage))], 'trap_ID': [trap_id for _ in np.arange(len(lineage))], 'trace': [trace for _ in np.arange(len(lineage))],
    #             'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
    #         }
    #         to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables_new})
    #         to_add = pd.DataFrame(to_add)
    #         means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    #
    #
    #     # for dataset in ['SL', 'NL']:
    #     #     for trap_id in np.unique(info[(info['dataset'] == dataset)]['trap_ID']):
    #     #         for trace in ['A', 'B']:
    #     #
    #     #             # the values of the lineage we get from physical units
    #     #             lineage = info[(info['trap_ID'] == trap_id) & (info['dataset'] == dataset) & (info['trace'] == trace)].copy()
    #     #
    #     #             # add its time-average
    #     #             to_add = {
    #     #                 'dataset': [dataset for _ in np.arange(len(lineage))], 'trap_ID': [trap_id for _ in np.arange(len(lineage))], 'trace': [trace for _ in np.arange(len(lineage))],
    #     #                 'max_gen': [len(lineage) for _ in np.arange(len(lineage))], 'generation': np.arange(len(lineage))
    #     #             }
    #     #             to_add.update({param: [np.mean(lineage[param]) for _ in np.arange(len(lineage))] for param in phenotypic_variables_new})
    #     #             to_add = pd.DataFrame(to_add)
    #     #             means_df = means_df.append(to_add, ignore_index=True).reset_index(drop=True)
    #
    #     assert len(info) == len(means_df)
    
    return means_df


""" Creates a dataframe with the ergodicity breaking parameter of each phenotypic variable """


def ergodicity_breaking_parameter(df, phenotypic_variables_new, kind, n_boots=0):  # MM
    # Initialize where we will put the bootstrapped ergodicity breaking variable
    eb_df = pd.DataFrame(columns=['variable', 'kind', 'value'])
    
    # get the dataframes where in each entry, instead of the generation specific value, there is the time-average
    time_averages = get_time_averages_df(df, phenotypic_variables_new)  # MM
    
    # to normalize the different phenotypic_variables_new
    pop_var = df.var()
    
    # bootstrap this ergodicity breaking parameter
    if n_boots != 0:
        for _ in np.arange(n_boots):
            # Bootstrapping the indices has the lineage length inequality taken into account
            indices = time_averages.sample(frac=1, replace=True).index
            
            # Get the variance of the time-averages for both kinds of lineages
            variance = time_averages[phenotypic_variables_new].loc[indices].var() / pop_var
            
            # add them both to the dataframe where we save it all
            for param in phenotypic_variables_new:
                eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    else:
        # Get the variance of the time-averages for both kinds of lineages
        variance = time_averages[phenotypic_variables_new].var() / pop_var
    
        # add them both to the dataframe where we save it all
        for param in phenotypic_variables_new:
            eb_df = eb_df.append({'variable': param, 'kind': kind, 'value': variance[param]}, ignore_index=True)
    
    return [time_averages, eb_df]


""" Creates this dataframe for all kinds of data we have """


def main(args):
    
    """
    Creates the dataframes with the
    :param args: physical units, save folder + names of all the returning dfs
    :return: the population sampled df, time_averages, ergodicity breaking parameter df, kl divergences df
    """
    
    for blah in ['A', 'B']:
        # import the labeled measured bacteria in physical units
        info = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
        
        #HERE WE GO
        info = info[info['trace'] == blah].copy()
    
        # Add what Naama wanted
        size_fraction = np.array([])
    
        for lin_id in info.lineage_ID.unique():
            lineage = info[info['lineage_ID'] == lin_id].copy().sort_values('generation')
            fr = lineage['length_birth'].values[1:] / lineage['length_birth'].values[:-1]
            fr = np.append(fr, np.array([np.nan]))
            size_fraction = np.append(size_fraction, fr)
    
        info['size_fraction'] = pd.Series(size_fraction.flatten())
        
        print('Population physical_units lineages')
        # Create lineages sampled from a population distribution
        shuffled = shuffle_info(info, args.MM)
        # shuffled.to_csv('{}/{}'.format(args.save_folder, args.population_sampled), index=False)
        
        print('ergodicity breaking and time-averages')
        # get the bootstrapped EB variable for both kinds of lineages
        time_averages_trace, eb_df = ergodicity_breaking_parameter(info, phenotypic_variables_new, kind='Trace')
        _, eb_df_pop = ergodicity_breaking_parameter(shuffled, phenotypic_variables_new, kind='Artificial')
        eb_df = eb_df.append(eb_df_pop, ignore_index=True).reset_index(drop=True)
        
        # save it to the right folder
        # time_averages_trace.to_csv('{}/{}'.format(args.save_folder, args.ta), index=False)
        eb_df.to_csv('{}/{}_{}'.format(args.save_folder, args.ebp, blah), index=False)
        
        # print('kl_divergences')
        # # Put in the kl divergences for each variable for each type of lineage
        # kl_df = kl_divergence(info, phenotypic_variables_new, 'Trace')
        # kl_df = kl_df.append(kl_divergence(shuffled, phenotypic_variables_new, 'Artificial'), ignore_index=True).reset_index(drop=True)
        #
        # # save the kl_df dataframe
        # kl_df.to_csv('{}/{}'.format(args.save_folder, args.kld), index=False)
    
    
if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    # # Do all the Mother Machine data
    # for data_origin in dataset_names:
    #
    #     parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
    #     parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    #     parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
    #                         required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    #     parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
    #                         required=False, default='physical_units_with_outliers.csv')
    #     parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
    #                         required=False, default='population_lineages_with_outliers.csv')
    #     parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
    #                         required=False, default='time_averages_with_outliers.csv')
    #     parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
    #                         required=False, default='ergodicity_breaking_parameter_with_outliers.csv')
    #     parser.add_argument('-kld', '--kld', metavar='', type=str,
    #                         help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
    #                         required=False, default='kullback_leibler_divergences_with_outliers.csv')
    #     parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
    #     args = parser.parse_args()
    #
    #     main(args)
    #
    #     print('*' * 200)
        
    # How long did it take to do the mother machine?
    # mm_time = time.time() - first_time

    data_origin = 'SM'

    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units_with_outliers.csv')
    parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                        required=False, default='population_lineages_with_outliers.csv')
    parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                        required=False, default='time_averages_with_outliers.csv')
    parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                        required=False, default='ergodicity_breaking_parameter_with_outliers.csv')
    parser.add_argument('-kld', '--kld', metavar='', type=str,
                        help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                        required=False, default='kullback_leibler_divergences_with_outliers.csv')
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    args = parser.parse_args()
    
    main(args)

    # # How long did it take to do the mother machine?
    # sm_time = time.time() - (mm_time + first_time)
    #
    # print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
