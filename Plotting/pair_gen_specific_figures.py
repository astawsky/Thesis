#!/usr/bin/env bash

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from string import ascii_uppercase as uppercase_letters
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets, seaborn_preamble


def plot_panel(var_a, var_b, count1, pair_correlation, ax):
    cmap = sns.color_palette('tab10')[5:8]
    cmap_dark = sns.color_palette('pastel')[5:8]
    
    # The constant
    df = pair_correlation[(pair_correlation['kind'] == 'time_averages') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)].copy()

    # if count1 == 0:
    #     sns.catplot(data=pair_correlation[(pair_correlation['kind'] == 'trace_centered') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)], x='generation',
    #                 y='correlation', hue='dataset', palette=cmap, edgecolor='black', ax=ax, legend=True, kind='bar')
    # else:
    #     sns.catplot(data=pair_correlation[(pair_correlation['kind'] == 'trace_centered') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)], x='generation',
    #                 y='correlation', hue='dataset', palette=cmap, edgecolor='black', ax=ax, legend=False, kind='bar')

    if count1 == 0:
        sns.barplot(data=pair_correlation[(pair_correlation['kind'] == 'trace_centered') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)], x='generation',
                    y='correlation', hue='dataset', palette=cmap, edgecolor='black', ax=ax)
    else:
        sns.barplot(data=pair_correlation[(pair_correlation['kind'] == 'trace_centered') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)], x='generation',
                    y='correlation', hue='dataset', palette=cmap, edgecolor='black', ax=ax)
        ax.get_legend().remove()

    # sns.barplot(data=pair_correlation[(pair_correlation['kind'] == 'trace_centered') & (pair_correlation['variable_a'] == var_a) & (pair_correlation['variable_b'] == var_b)], x='generation',
    #             y='correlation', hue='dataset', palette=cmap, edgecolor='black', ax=ax, legend_out=True)
    # ax.get_legend().remove()
    
    # Plot the constant impact of the time-average in the
    for count, hue in enumerate(datasets):
        
        # first line is for a makeshift border, otherwise it is difficult to see the horizontal lines with the same colors as the barplots
        # ax.axhline(np.unique(df[df['dataset'] == hue]['correlation'])+.01, ls='--', color='black', linewidth=2)
        # ax.axhline(np.unique(df[df['dataset'] == hue]['correlation'])-.01, ls='--', color='black', linewidth=2)
        if count1 == 1:
            ax.axhline(np.unique(df[df['dataset'] == hue]['correlation']), ls='--', color=cmap_dark[count], linewidth=3, label=hue)
        else:
            ax.axhline(np.unique(df[df['dataset'] == hue]['correlation']), ls='--', color=cmap_dark[count], linewidth=3)
    
    
def plot_main_figure(args, pair_correlation):
    seaborn_preamble()
    fig, axes = plt.subplots(2, 2, tight_layout=True, figsize=[7, 7], sharey=True, sharex=True)
    axes = axes.flatten()
    for count1, var_a in enumerate(['fold_growth', 'generationtime', 'length_birth', 'growth_rate']):
        var_b = var_a
        ax = axes[count1]
        plot_panel(var_a, var_b, count1, pair_correlation, ax)
        ax.set_xlabel('')
        ax.set_ylabel(symbols['trace_centered'][var_a], fontsize='xx-large')

        # Let's put the titles
        ax.set_title(uppercase_letters[count1], x=-0.07, fontsize='xx-large')
        
        # the legend
        if count1 == 0:
            ax.legend(title='')
        elif count1 == 1:
            h, l = ax.get_legend_handles_labels()
            ax.legend(h[:3], l[:3], title='')
        # if count1 in [0, 1]:
        #     ax.legend(title='')
    
    plt.savefig('{}/gen_specific_main_figure.png'.format(args.figs_location), dpi=300)
    # plt.show()
    plt.close()
    

# def plot_three_in_one(args):
#     repeats = []
#     for var_a in phenotypic_variables:
#         for var_b in phenotypic_variables:
#             if var_b not in repeats:
#
#                 seaborn_preamble()
#                 fig, ax = plt.subplots(1, 1, tight_layout=True, figsize=[7, 7])
#
#                 plot_panel(var_a, var_b, ax)
#
#                 plt.legend()
#                 plt.ylabel('')
#                 plt.xlabel('')
#                 plt.title('{} and {}'.format(symbols['physical_units'][var_a], symbols['physical_units'][var_b]))
#                 if var_a == var_b:
#                     plt.savefig('{}/Diagonal gen-specific correlations/{}'.format(args.figs_location, var_a), dpi=300)
#                 else:
#                     plt.savefig('{}/Non-Diagonal gen-specific correlations/{}, {}'.format(args.figs_location, var_a, var_b), dpi=300)
#                 # plt.show()
#                 plt.close()
#         repeats.append(var_a)
        
        
def main(args):
    create_folder(args.figs_location + '/Diagonal gen-specific correlations')
    create_folder(args.figs_location + '/Non-Diagonal gen-specific correlations')
    
    # save it to the Data folder
    pair_correlation = pd.read_csv('{}/over_lineages.csv'.format(args.save_folder))
    
    plot_main_figure(args, pair_correlation)
