#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, shuffle_info, dataset_names, seaborn_preamble, symbols, create_folder
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


# """ adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """
#
#
# def expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df, label):
#     import matplotlib.pyplot as plt
#
#     # create the info dataframe with values being the cumulative mean of the lineage
#     for dataset in np.unique(df['dataset']):
#
#         for trap_id in np.unique(df[df['dataset'] == dataset]['trap_ID']):
#             # Go through both traces in an SM trap
#             for trace_id in ['A', 'B']:
#                 # specify the trace
#                 trace = df[(df['trap_ID'] == int(trap_id)) & (df['trace'] == trace_id) & (df['dataset'] == dataset)].copy()
#
#                 # add its time-average up until and including this generation
#                 to_add = pd.DataFrame.from_dict({
#                     'label': [label for _ in range(len(trace))],
#                     'trap_ID': [int(trap_id) for _ in range(len(trace))],
#                     'trace': [trace_id for _ in range(len(trace))],
#                     'generation': [generation + 1 for generation in trace['generation']]
#                 }).reset_index(drop=True)
#
#                 # plt.axhline(0, color='black')
#                 # plt.plot(trace.sort_values('generation')['generationtime'].expanding().mean().values, label='expanding mean')
#                 # plt.plot(trace.sort_values('generation')['generationtime'].cumsum().values, label='cumsum')
#                 # plt.legend()
#                 # plt.show()
#                 # plt.close()
#                 # exit()
#
#                 to_add_expanding_mean = pd.concat([trace.sort_values('generation')[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
#                 to_add_cum_sum = pd.concat([trace.sort_values('generation')[phenotypic_variables].cumsum().reset_index(drop=True), to_add], axis=1)
#                 expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
#                 cum_sum_df = cum_sum_df.append(to_add_cum_sum, ignore_index=True).reset_index(drop=True)
#
#                 if to_add_expanding_mean.isnull().values.any():
#                     print('expanding mean')
#                     print(dataset, trap_id, trace)
#                     print(to_add_expanding_mean)
#                     input()
#                 if cum_sum_df.isnull().values.any():
#                     print('expanding mean')
#                     print(dataset, trap_id, trace)
#                     print(cum_sum_df)
#                     input()
#
#     assert not expanding_mean.isnull().values.any()
#     assert not cum_sum_df.isnull().values.any()
#
#     # Calculate the cv and var over all lineages in the dataset
#     for param in phenotypic_variables:
#         for generation in np.arange(1, df['generation'].max()+1):
#             time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)][param]
#             different_walks = cum_sum_df[(cum_sum_df['label'] == label) & (cum_sum_df['generation'] == generation)][param]
#
#             if time_averages.isnull().values.any():
#                 print('expanding mean')
#                 print(time_averages)
#                 input()
#             if different_walks.isnull().values.any():
#                 print('expanding mean')
#                 print(different_walks)
#                 input()
#
#             # print(time_averages)
#             # print(different_walks)
#             # print(time_averages.isnull().values.any())
#             # print(different_walks.isnull().values.any())
#             # exit()
#
#             # adding it to the dataframe of all coefficients of variations
#             variance_of_expanding_mean = variance_of_expanding_mean.append(
#                 {
#                     'label': label,
#                     'generation': generation,
#                     'cv': time_averages.std() / time_averages.mean(),
#                     'var': time_averages.var(),
#                     'param': param
#                 }, ignore_index=True
#             )
#             variance_of_cum_sum_df = variance_of_cum_sum_df.append(
#                 {
#                     'label': label,
#                     'generation': generation,
#                     'var': different_walks.var(),
#                     'param': param
#                 }, ignore_index=True
#             )
#
#     return [expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df]
#
#
# def main(args):
#     # # import/create the trace lineages
#     # physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
#     #
#     # physical_units = limit_lineage_length(physical_units, min_gens=50)
#     # trace_centered = trace_center_a_dataframe(physical_units)
#
#     # print(physical_units.isnull().any())
#     # print(trace_centered.isnull().any())
#     # print(len(physical_units), len(trace_centered))
#     # exit()
#
#     # import/create the trace lineages
#     physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
#
#     # import/create the trace-centered lineages
#     trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tc))
#
#     # import/create the population lineages
#     try:
#         population_sampled = pd.read_csv('{}/{}'.format(args.save_folder, args.population_sampled))
#     except:
#         print('creating population sampled')
#         population_sampled = shuffle_info(physical_units, mm=args.MM)
#         population_sampled.to_csv('{}/{}'.format(args.save_folder, args.population_sampled), index=False)
#
#     # import/create the shuffled generations lineages
#     try:
#         shuffled_generations = pd.read_csv('{}/{}'.format(args.save_folder, args.shuffled))
#     except:
#         print('creating shuffled generations')
#         shuffled_generations = shuffle_lineage_generations(physical_units)
#         shuffled_generations.to_csv('{}/{}'.format(args.save_folder, args.shuffled), index=False)
#
#     # The generation-shuffled trace-centered dataframe
#     try:
#         shuffled_tc = pd.read_csv('{}/{}'.format(args.save_folder, args.tc_shuffled))
#     except:
#         print('creating shuffled tc')
#         shuffled_tc = shuffle_lineage_generations(trace_centered)
#         shuffled_tc.to_csv('{}/{}'.format(args.save_folder, args.tc_shuffled), index=False)
#
#     # We keep the trap means here
#     expanding_mean = pd.DataFrame(columns=['label', 'trap_ID', 'trace', 'generation'] + phenotypic_variables)
#
#     # Keep the cv per lineage length here here
#     variance_of_expanding_mean = pd.DataFrame(columns=['label', 'generation', 'cv', 'var', 'param'])
#
#     # We keep the trap means here
#     cum_sum_df = pd.DataFrame(columns=['label', 'trap_ID', 'trace', 'generation'] + phenotypic_variables)
#
#     # Keep the cv per lineage length here here
#     variance_of_cum_sum_df = pd.DataFrame(columns=['label', 'generation', 'var', 'param'])
#
#     # Calculate the cv and TA per lineage length
#     for kind, df in zip(['Shuffled', 'Trace', 'Population', 'Trace-Centered', 'Shuffled TC'], [shuffled_generations, physical_units, population_sampled, trace_centered, shuffled_tc]):
#         # kind = 'Trace-Centered'
#         # df = trace_centered
#         print(kind)
#         expanding_mean, variance_of_expanding_mean, cum_sum_df, variance_of_cum_sum_df = expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, variance_of_expanding_mean,
#                                                                                                                              cum_sum_df,
#                                                                                                                              variance_of_cum_sum_df, kind)
#
#     # Save the csv file
#     expanding_mean.to_csv('{}/{}'.format(args.save_folder, args.cta), index=False)
#
#     # Save the csv file
#     variance_of_expanding_mean.to_csv('{}/{}'.format(args.save_folder, args.vcta), index=False)
#
#     # Save the csv file
#     cum_sum_df.to_csv('{}/{}'.format(args.save_folder, args.cum_sum), index=False)
#
#     # Save the csv file
#     variance_of_cum_sum_df.to_csv('{}/{}'.format(args.save_folder, args.v_cum_sum), index=False)


# """ adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """
#
#
# def expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, label, out_df, first_generation, final_generation):
#     pooled_ensemble_total = pd.DataFrame(columns=df.columns)
#     print(df.columns)
#     # exit()
#
#     # create the info dataframe with values being the cumulative mean of the lineage
#     for lin_id in df.lineage_ID.unique():
#
#         # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
#         lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
#
#         # add its time-average up until and including the total number of generations without NaNs
#         to_add = pd.DataFrame.from_dict({
#             'label': [label for _ in range(len(lineage))],
#             'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
#             'generation': [generation for generation in np.arange(len(lineage), dtype=int)]
#         }).reset_index(drop=True)
#
#         # Get the expanding mean and the cumulative sum, make sure they are robust to the first NaN in division ratio
#         to_add_expanding_mean = pd.concat([lineage[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
#         expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
#
#         if to_add_expanding_mean.isnull().values.any():
#             raise IOError('Got NaNs in the expanding mean, meaning something is wrong')
#
#         lineage['generation'] = np.arange(len(lineage), dtype=int)
#
#         pooled_ensemble_total = pooled_ensemble_total.append(lineage, ignore_index=True)
#
#         if len(to_add_expanding_mean) != len(lineage):
#
#             print(lineage)
#             print(to_add_expanding_mean)
#             raise IOError('Wronfg lenghts')
#
#         if not np.array_equal(to_add_expanding_mean['generation'].values, lineage['generation'].values):
#             print(lin_id)
#             print(type(lineage['generation'].iloc[0]))
#             print(type(to_add_expanding_mean['generation'].iloc[0]))
#             raise IOError('Wronfg values of generation')
#
#     # This is important
#     assert not expanding_mean.isnull().values.any()
#
#     # This is important
#     assert not pooled_ensemble_total.isnull().values.any()
#
#     # Calculate the cv and var over all lineages in the dataset
#     for generation in np.arange(first_generation, final_generation + 1):
#         print(generation)
#         generation = 30
#
#         # time_averages1 = [
#         #     df[(df['lineage_ID'] == lin_id) & (df['generation'] <= generation)].sort_values('generation').dropna(axis=0)[phenotypic_variables].mean()
#         #     for lin_id in df.lineage_ID.unique() if len(df[(df['lineage_ID'] == lin_id)].dropna(axis=0)) >= generation
#         # ]
#
#         # Get the time-averages of the lineages that have the generation we are looking for
#         time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)].copy()  # Excludes all lineages that do not reach this generation
#
#         # print(time_averages, len(time_averages))
#         # print(time_averages1, len(time_averages1))
#         # exit()
#
#         # Get the array of which lineages are long enough, ie. have the amount of cycles we need
#         long_enough_lineages = time_averages['lineage_ID'].unique()
#
#         # Define the pooled ensemble of all the lineages that have at least this generation
#         pooled_ensemble = pooled_ensemble_total[(pooled_ensemble_total['generation'] <= generation) & (pooled_ensemble_total['lineage_ID'].isin(long_enough_lineages))].copy()
#         pooled_ensemble1 = df.dropna(axis=0)[(df['generation'] <= generation) & (df['lineage_ID'].isin(long_enough_lineages))].copy()
#
#         # print(len(pooled_ensemble), len(pooled_ensemble1))
#         # print(pooled_ensemble)
#         # print(pooled_ensemble['dataset'].unique())
#         # print(pooled_ensemble1['dataset'].unique())
#         # exit()
#
#         # That way we have a robust variance across lineages
#         if len(long_enough_lineages) > 1:
#
#             # for each parameter combination decompose the pooled variance
#             repeat = []
#             for param1 in phenotypic_variables:
#                 for param2 in phenotypic_variables:
#                     if param2 not in repeat:
#
#                         print(param1, param2)
#
#                         # The mean is a linear operator
#                         if pooled_ensemble[param2].mean() != time_averages[param2].mean():
#                             print(pooled_ensemble[param2].mean(), time_averages[param2].mean())
#
#                         print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
#                         exit()
#
#                         # We must have the same amount of
#                         if not (generation * len(time_averages[param1])) == (generation * len(time_averages[param2])) == len(pooled_ensemble[param1]) == len(pooled_ensemble[param2]):
#                             print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
#                             exit()
#
#                         print('passed!')
#                         exit()
#
#                         # Get the variances of each type of series
#                         ta_cov = (generation * ((time_averages[param1].copy() - time_averages.mean()[param1].copy()) * (time_averages[param2].copy() - time_averages.mean()[param2].copy()))).sum() / (len(pooled_ensemble) - 1)
#                         pool_var = ((pooled_ensemble[param1].copy() - pooled_ensemble.mean()[param1].copy()) * (pooled_ensemble[param2].copy() - pooled_ensemble.mean()[param2].copy())).sum() / (len(pooled_ensemble) - 1)
#
#                         # Calculate the variance of the time-averages
#                         gamma_ta_cov = ta_cov / pool_var
#
#                         print(gamma_ta_cov)
#
#                         # Add it to the dataframe to output
#                         out_df = out_df.append({
#                             'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
#                         }, ignore_index=True)
#
#                 repeat.append(param1)
#
#         exit()
#
#     print(out_df)
#
#     # Important
#     if out_df.isnull().values.any():
#         print('ERROR!')
#         print(out_df)
#         exit()
#
#     return out_df


""" adds to one inputted dataframe the cumulative means, and to the the other inputed dataframe the variance and cv of the TAs of every cycle parameter per lineage length """


def expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, label, out_df, first_generation, final_generation):
    # For the pooled ensemble
    pooled_ensemble_total = pd.DataFrame(columns=df.columns)
    
    # create the info dataframe with values being the cumulative mean of the lineage
    for lin_id in df.lineage_ID.unique():
        
        # specify the trace, we drop all the NaN rows/cycles because we need the same number of samples for all variables in order to get the covariance
        lineage = df[(df['lineage_ID'] == lin_id)].sort_values('generation').dropna(axis=0).copy().reset_index(drop=True)
        
        # add its time-average up until and including the total number of generations without NaNs
        to_add = pd.DataFrame.from_dict({
            'label': [label for _ in range(len(lineage))],
            'lineage_ID': [int(lin_id) for _ in range(len(lineage))],
            'generation': [generation for generation in np.arange(len(lineage), dtype=int)]
        }).reset_index(drop=True)
        
        # Get the expanding mean and the cumulative sum, make sure they are robust to the first NaN in division ratio
        to_add_expanding_mean = pd.concat([lineage[phenotypic_variables].expanding().mean().reset_index(drop=True), to_add], axis=1)
        expanding_mean = expanding_mean.append(to_add_expanding_mean, ignore_index=True).reset_index(drop=True)
        
        if to_add_expanding_mean.isnull().values.any():
            raise IOError('Got NaNs in the expanding mean, meaning something is wrong')
        
        lineage['generation'] = np.arange(len(lineage), dtype=int)
        
        pooled_ensemble_total = pooled_ensemble_total.append(lineage, ignore_index=True)
        
        if len(to_add_expanding_mean) != len(lineage):
            
            print(lineage)
            print(to_add_expanding_mean)
            raise IOError('Wronfg lenghts')
        
        if not np.array_equal(to_add_expanding_mean['generation'].values, lineage['generation'].values):
            print(lin_id)
            print(type(lineage['generation'].iloc[0]))
            print(type(to_add_expanding_mean['generation'].iloc[0]))
            raise IOError('Wronfg values of generation')
    
    # This is important
    assert not expanding_mean.isnull().values.any()
    
    # This is important
    assert not pooled_ensemble_total.isnull().values.any()
    
    # Calculate the cv and var over all lineages in the dataset
    for generation in np.arange(first_generation, final_generation + 1):
        print(generation)
        
        # Get the time-averages of the lineages that have the generation we are looking for
        time_averages = expanding_mean[(expanding_mean['label'] == label) & (expanding_mean['generation'] == generation)].copy()  # Excludes all lineages that do not reach this generation
        
        # print(time_averages[['generationtime', 'growth_rate']])
    
        # Get the array of which lineages are long enough, ie. have the amount of cycles we need
        long_enough_lineages = time_averages['lineage_ID'].unique()
    
        # Define the pooled ensemble of all the lineages that have at least this generation
        pooled_ensemble = pooled_ensemble_total[(pooled_ensemble_total['generation'] <= generation) & (pooled_ensemble_total['lineage_ID'].isin(long_enough_lineages))].copy()
        
        # That way we have a robust variance across lineages
        if len(long_enough_lineages) > 1:
            
            # for each parameter combination decompose the pooled variance
            repeat = []
            for param1 in phenotypic_variables:
                for param2 in phenotypic_variables:
                    if param2 not in repeat:
                        
                        # The mean is a linear operator
                        if np.abs(pooled_ensemble[param2].mean() - time_averages[param2].mean()) > 0.000001:
                            print(pooled_ensemble[param2].mean(), time_averages[param2].mean())
                            print(pooled_ensemble[param2].mean() - time_averages[param2].mean())
                            raise IOError('Means are not the same in a big way')
                        
                        # We must have the same amount of
                        if not ((generation + 1) * len(time_averages[param1])) == ((generation + 1) * len(time_averages[param2])) == len(pooled_ensemble[param1]) == len(pooled_ensemble[param2]):
                            print((generation * len(time_averages[param1])), (generation * len(time_averages[param2])), len(pooled_ensemble[param1]), len(pooled_ensemble[param2]))
                            raise IOError('Sizes are not the same in a big way')
        
                        # Get the variances of each type of series
                        ta_cov = ((generation + 1) * ((time_averages[param1].copy() - time_averages.mean()[param1].copy()) * (time_averages[param2].copy() - time_averages.mean()[param2].copy()))).sum() / (len(pooled_ensemble) - 1)
                        pool_var = pooled_ensemble[param1].std() * pooled_ensemble[param2].std()
        
                        # Calculate the variance of the time-averages
                        gamma_ta_cov = ta_cov / pool_var
                
                        # Add it to the dataframe to output
                        out_df = out_df.append({
                            'param1': param1, 'param2': param2, 'n': len(long_enough_lineages), 'generation': generation, 'gamma_ta': gamma_ta_cov, 'label': label
                        }, ignore_index=True)
                    
                repeat.append(param1)
        
    print(out_df)

    # Important
    if out_df.isnull().values.any():
        print('ERROR!')
        print(out_df)
        exit()
    
    return out_df


def main(args, first_generation, final_generation):
    # import/create the trace lineages
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    # if args.data_origin == 'SM':
    #     # import/create the trace lineages
    #     physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).dropna()
    # else:
    #     # import/create the trace lineages
    #     physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))

    population_sampled = shuffle_info(physical_units, mm=args.MM)
    
    # # import/create the population lineages
    # try:
    #     if args.data_origin == 'SM':
    #         population_sampled = pd.read_csv('{}/{}'.format(args.save_folder, args.population_sampled)).dropna()
    #     else:
    #         population_sampled = pd.read_csv('{}/{}'.format(args.save_folder, args.population_sampled))
    # except:
    #     print('creating population sampled')
    #     population_sampled = shuffle_info(physical_units, mm=args.MM)
    #     population_sampled.to_csv('{}/{}'.format(args.save_folder, args.population_sampled), index=False)
    #
    #     if args.data_origin == 'SM':
    #         population_sampled = population_sampled.dropna()
    
    # We keep the trap means here
    expanding_mean = pd.DataFrame(columns=['label', 'lineage_ID', 'generation'] + phenotypic_variables)

    # Where we keep the gammas
    out_df = pd.DataFrame(columns=['param1', 'param2', 'n', 'generation', 'gamma_ta', 'label'])
    
    # Calculate the cv and TA per lineage length
    for kind, df in zip(['Trace', 'Artificial'], [physical_units, population_sampled]):
        
        print(kind)
        out_df = expanding_mean_cumsum_and_variances(df, phenotypic_variables, expanding_mean, kind, out_df, first_generation, final_generation)
        
    out_df.to_csv('{}/{}'.format(args.save_folder, 'gamma_ta_corrs_per_gen.csv'), index=False)
    
    
def plot(args):
    import seaborn as sns
    
    gamma_ta_corrs_per_gen = pd.read_csv('{}/{}'.format(args.save_folder, 'gamma_ta_corrs_per_gen.csv')).sort_values('generation')

    seaborn_preamble()
    plt.plot(gamma_ta_corrs_per_gen['generation'].values, gamma_ta_corrs_per_gen['n'].values)
    plt.title('lineages per generation')
    # plt.show()
    plt.tight_layout()
    plt.savefig('{}/{}/lineages per generation.png'.format(args.figs_location, args.sc_decomposition_per_gen), dpi=300)
    plt.close()
    
    repeat = []
    for param1 in gamma_ta_corrs_per_gen.param1.unique():
        for param2 in gamma_ta_corrs_per_gen.param2.unique():
            if param2 not in repeat:
                print(param1, param2)
                
                # to_plot['index'] = list(to_plot.index)
                # to_plot = to_plot[(to_plot['index'].isin(to_plot.index.values[:2]))]
                # print(to_plot)
                # print(to_plot['generation'].values, to_plot['gamma_ta'].values, sep='\n')
                # sns.pointplot(to_plot['generation'].values, to_plot['gamma_ta'].values)
                # plt.axhline(0)
                for label in ['Trace', 'Artificial']:
                    
                    # param1 = 'growth_rate'
                    # param2 = 'length_birth'
                    
                    # What are the values of Gamma TA per gen
                    to_plot = gamma_ta_corrs_per_gen[(gamma_ta_corrs_per_gen['param1'] == param1) & (gamma_ta_corrs_per_gen['param2'] == param2) & (gamma_ta_corrs_per_gen['label'] == label)].sort_values('generation').copy()
                    
                    # Wrong ordering of the variables for some reason
                    if len(to_plot) == 0:
                        to_plot = gamma_ta_corrs_per_gen[(gamma_ta_corrs_per_gen['param1'] == param2) & (gamma_ta_corrs_per_gen['param2'] == param1) & (gamma_ta_corrs_per_gen['label'] == label)].sort_values('generation').copy()
                        
                    if param1 == param2:
                        slope, intercept = linregress(np.log(to_plot.generation.values), np.log(to_plot.gamma_ta.values))[:2]
                        sns.scatterplot(data=to_plot, x='generation', y='gamma_ta', label=label + r': $\Gamma = '+str(np.exp(intercept))[:4]+' \cdot n^{'+str(slope)[:4]+'}$')
                        plt.plot(to_plot.generation.unique(), [np.exp(intercept)*(l**slope) for l in to_plot.generation.unique()])
                        plt.yscale('log')
                        plt.xscale('log')
                    else:
                        sns.scatterplot(data=to_plot, x='generation', y='gamma_ta', label=label)
                plt.legend()
                if param1 == param2:
                    plt.title(r'$\Gamma($'+symbols['with_n_ta'][param1]+r'$)$')
                    plt.savefig('{}/{}/{} EB_n'.format(args.figs_location, args.sc_decomposition_per_gen, param1), dpi=300)
                else:
                    plt.title(r'$\Gamma($'+symbols['with_n_ta'][param1] + r'$, \, $' + symbols['with_n_ta'][param2] + r'$)$')
                    plt.savefig('{}/{}/{} {} EB_n'.format(args.figs_location, args.sc_decomposition_per_gen, param1, param2), dpi=300)
                # plt.show()
                plt.close()
        
        repeat.append(param1)
    
    
if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # The final_generations for the different datasets
    fg = {
        'SM': 45,
        'lambda_LB': 60,
        'Maryam_LongTraces': 45,
        'MG1655_inLB_LongTraces': 200,
        'LAC_M9': np.nan
    }
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in dataset_names:
        print(data_origin)
        if data_origin == 'lambda_LB':
            continue

        # # This is because this dataset does not have enough generations for the analysis to be relevant
        # if data_origin == 'MG1655_inLB_LongTraces':  # data_origin == 'LAC_M9' or data_origin == 'lambda_LB' or
        #     continue

        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                            required=False, default='population_lineages.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                            required=False, default='trace_centered.csv')
        parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        parser.add_argument('-sc_decomposition_per_gen', '--sc_decomposition_per_gen', metavar='', type=str, help='The name of the single cell decomposition figures.',
                            required=False, default='sc_decomposition_per_gen')
        parser.add_argument('-lineages_per_gen', '--lineages_per_gen', metavar='', type=str, help='The name of the lineages_per_gen figures.',
                            required=False, default='lineages_per_gen')
        args = parser.parse_args()

        create_folder('{}/{}'.format(args.figs_location, args.sc_decomposition_per_gen))

        main(args, first_generation=5, final_generation=fg[data_origin])

        plot(args)

        print('*' * 200)

    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time

    data_origin = 'SM'

    print(data_origin)
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                        required=False, default='population_lineages.csv')
    parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                        required=False, default='time_averages.csv')
    parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                        required=False, default='trace_centered.csv')
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    parser.add_argument('-sc_decomposition_per_gen', '--sc_decomposition_per_gen', metavar='', type=str, help='The name of the single cell decomposition figures.',
                        required=False, default='sc_decomposition_per_gen')
    parser.add_argument('-lineages_per_gen', '--lineages_per_gen', metavar='', type=str, help='The name of the lineages_per_gen figures.',
                        required=False, default='lineages_per_gen')
    args = parser.parse_args()

    create_folder('{}/{}'.format(args.figs_location, args.sc_decomposition_per_gen))

    main(args, first_generation=5, final_generation=fg[data_origin])

    plot(args)

    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)

    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
