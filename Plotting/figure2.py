#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import symbols, create_folder, seaborn_preamble, dataset_names
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
            

# """ Create the dataframe that contains the regression phenotypic_variables """
#
#
# def slope_comparisons(args, variation_of_expanding_mean):
#
#     # This is where we will keep the regression phenotypic_variables
#     slope_df = pd.DataFrame(columns=['Trace', 'Shuffled', 'Population', 'Trace-Centered'])
#     intercept_df = slope_df.copy()
#
#     # looping over the phenotypic_variables in highlights and their respective axes
#     for count, param in enumerate(phenotypic_variables):
#         to_add_intercept = {}
#         to_add_slope = {}
#
#         # Plot the four different kinds of lineages
#         for label, marker in zip(slope_df.columns, ['o', 'x', '^', '+']):
#
#             # the variance and generations of the TAs of this type of lineage
#             dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
#
#             # So each kind of lineage starts at the 1 and keeps their slope
#             dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
#
#             # We get the best fit line
#             lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
#
#             # Add this regression phenotypic_variables to the dataframe
#             to_add_slope.update({label: lin_reg.coef_[0]})
#             to_add_intercept.update({label: lin_reg.intercept_[0]})
#
#         # Add it to the dataframe to make a table or heatmap
#         intercept_df = intercept_df.append(pd.DataFrame(to_add_intercept, index=[symbols['time_averages'][param]]))
#         slope_df = slope_df.append(pd.DataFrame(to_add_slope, index=[symbols['time_averages'][param]]))
#
#     # set a style on seaborn module
#     seaborn_preamble()
#
#     fig, ax = plt.subplots(figsize=[7, 3.5])
#
#     sns.heatmap(slope_df.T, annot=True, square=True, yticklabels=True)
#     # plt.title('Slopes')
#     ax.set_yticklabels(['Trace', 'Shuffled', 'Population', 'Trace-Centered'], rotation=0)
#     plt.tight_layout()
#     # plt.show()
#     plt.savefig('{}/Slopes heatmap.png'.format(args.figs_location), dpi=300)
#     plt.close()
#
#     fig, ax = plt.subplots(figsize=[7, 3.5])
#
#     sns.heatmap(intercept_df.T, annot=True, square=True, yticklabels=True)
#     # plt.title('Intercepts')
#     ax.set_yticklabels(['Trace', 'Shuffled', 'Population', 'Trace-Centered'], rotation=0)
#     plt.tight_layout()
#     # plt.show()
#     plt.savefig('{}/Intercepts heatmap.png'.format(args.figs_location), dpi=300)
#     plt.close()
#
#
# """ This creates the big figure 2 from the paper. Includes the individual figures code where  """
#
#
# def variance_timeaverage_per_generation(args, variation_of_expanding_mean):
#
#     # The cycle phenotypic_variables we will look at
#     highlights = ['length_birth', 'growth_rate', 'generationtime']
#
#     # set a style on seaborn module
#     seaborn_preamble()
#
#     # create the figure and construct the layout of the figure
#     fig, axes = plt.subplots(nrows=1, ncols=len(highlights), sharey=True, figsize=[7, 3])
#
#     # looping over the phenotypic_variables in highlights and their respective axes
#     for count, param in enumerate(highlights):
#         # Specify the axis
#         ax = axes[count]
#
#         # Plot the four different kinds of lineages
#         for label, marker in zip(['Trace', 'Population', 'Trace-Centered', 'Shuffled'], ['o', 'x', '^', '+']):
#
#             # the variance and generations of the TAs of this type of lineage
#             dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
#
#             # So each kind of lineage starts at the 1 and keeps their slope
#             dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
#
#             # # We get the best fit line
#             # lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
#             # ax.plot(dist['generation'], np.exp(lin_reg.predict(np.array(np.log(dist['generation'])).reshape(-1, 1)).flatten()), ls='--', label='')
#
#             # plot the plots
#             ax.scatter(dist['generation'], dist['var'], marker=marker, alpha=.5)
#
#             # reg_order = PowerRegression(dist['generation'], dist['var'], bs=500)
#             # ax.scatter(dist['generation'], dist['var'], label=reg_order.param_labels['a'], marker=marker, alpha=.5)
#             # ax.plot(dist['generation'], reg_order.prediction, ls='--', label='')
#
#         # Add the grid
#         ax.grid(True)
#         # ax.legend()
#         ax.set(xscale="log", yscale="log")
#         # ax.tick_params(direction='out')
#         ax.set_xlabel(r'$N$')
#         ax.set_ylabel(r'$Var(${}$)$'.format(symbols['time_averages'][param], symbols['physical_units'][param]))
#         # ax.set_ylim([low_bound, high_bound])
#         # ax.set_xlim(left=0)
#
#         ax.set_title(uppercase_letters[count], x=-.15, fontsize='xx-large')
#
#     # Save the big figure 2.png
#     fig.tight_layout()
#     # plt.show()
#     plt.savefig('{}/Figure 2.png'.format(args.figs_location), dpi=300)
#     plt.close()
#
#
# """ This creates the rest of the figures left out of fig 2 in the main text to be in the appendix """
#
#
# def variance_timeaverage_per_generation_rest_of_figures(args, variation_of_expanding_mean):
#
#     # The cycle phenotypic_variables we will look at
#     highlights = ['fold_growth', 'length_final', 'division_ratio', 'added_length']
#
#     # set a style on seaborn module
#     seaborn_preamble()
#
#     # create the figure and construct the layout of the figure
#     fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=[7, 7])
#     axes = axes.flatten()
#
#     # looping over the phenotypic_variables in highlights and their respective axes
#     for count, param in enumerate(highlights):
#         # Specify the axis
#         ax = axes[count]
#
#         # Plot the four different kinds of lineages
#         for label, marker in zip(['Trace', 'Population', 'Trace-Centered', 'Shuffled'], ['o', 'x', '^', '+']):
#
#             # the variance and generations of the TAs of this type of lineage
#             dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
#
#             # So each kind of lineage starts at the 1 and keeps their slope
#             dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
#
#             # # We get the best fit line
#             # lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
#             # ax.plot(dist['generation'], np.exp(lin_reg.predict(np.array(np.log(dist['generation'])).reshape(-1, 1)).flatten()), ls='--', label='')
#
#             # plot the plots
#             ax.scatter(dist['generation'], dist['var'], marker=marker, alpha=.5)
#
#             # reg_order = PowerRegression(dist['generation'], dist['var'], bs=500)
#             # ax.scatter(dist['generation'], dist['var'], label=reg_order.param_labels['a'], marker=marker, alpha=.5)
#             # ax.plot(dist['generation'], reg_order.prediction, ls='--', label='')
#
#         # Add the grid
#         ax.grid(True)
#         # ax.legend()
#         ax.set(xscale="log", yscale="log")
#         # ax.tick_params(direction='out')
#         ax.set_xlabel(r'$N$')
#         ax.set_ylabel(r'$Var(${}$)$'.format(symbols['time_averages'][param]))
#         # ax.set_ylim([low_bound, high_bound])
#         # ax.set_xlim(left=0)
#
#         ax.set_title(uppercase_letters[count], x=-.15, fontsize='xx-large')
#
#     # Save the big figure 2.png
#     fig.tight_layout()
#     # plt.show()
#     plt.savefig('{}/Rest of Figure 2.png'.format(args.figs_location), dpi=300)
#     plt.close()


def main(args):
    import seaborn as sns
    
    gamma_ta_corrs_per_gen = pd.read_csv('{}/{}'.format(args.save_folder, 'gamma_ta_corrs_per_gen.csv'))

    seaborn_preamble()
    repeat = []
    for param1 in gamma_ta_corrs_per_gen.param1.unique():
        for param2 in gamma_ta_corrs_per_gen.param2.unique():
            if param2 not in repeat:
                # to_plot['index'] = list(to_plot.index)
                # to_plot = to_plot[(to_plot['index'].isin(to_plot.index.values[:2]))]
                # print(to_plot)
                # print(to_plot['generation'].values, to_plot['gamma_ta'].values, sep='\n')
                # sns.pointplot(to_plot['generation'].values, to_plot['gamma_ta'].values)
                # plt.axhline(0)
                for label in ['Trace', 'Artificial']:
                    to_plot = gamma_ta_corrs_per_gen[(gamma_ta_corrs_per_gen['param1'] == param1) & (gamma_ta_corrs_per_gen['param2'] == param2) & (gamma_ta_corrs_per_gen['label'] == label)].sort_values('generation').copy()
                    print(to_plot)
                    
                    slope, intercept = linregress(np.log(to_plot.generation.values), np.log(to_plot.gamma_ta.values))[:2]
                    sns.scatterplot(data=to_plot, x='generation', y='gamma_ta', label=label + r': $\Gamma = '+str(np.exp(intercept))[:4]+' \cdot n^{'+str(slope)[:4]+'}$')
                    plt.plot(to_plot.generation.unique(), [np.exp(intercept)*(l**slope) for l in to_plot.generation.unique()])
                plt.legend()
                plt.yscale('log')
                plt.xscale('log')
                plt.title(symbols['with_n'][param1])
                plt.savefig('{}/{} {} EB_n'.format(args.figs_location, param1, param2), dpi=300)
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
        'lambda_LB': 100,
        'Maryam_LongTraces': 70,
        'MG1655_inLB_LongTraces': 200,
        'LAC_M9': np.nan
    }
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in dataset_names:
        print(data_origin)
        
        # This is because this dataset does not have enough generations for the analysis to be relevant
        if data_origin == 'LAC_M9':
            continue
        
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
        args = parser.parse_args()
        
        create_folder('{}/{}'.format(args.figs_location, args.sc_decomposition_per_gen))
        
        main(args)
        exit()
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
    parser.add_argument('-sc_decomposition_per_gen', '--sc_decomposition_per_gen', metavar='', type=str, help='The name of the single cell decomposition per generation figures.',
                        required=False, default='sc_decomposition_per_gen')
    args = parser.parse_args()
    
    create_folder('{}/{}'.format(args.figs_location, args.sc_decomposition_per_gen))
    
    main(args)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))

