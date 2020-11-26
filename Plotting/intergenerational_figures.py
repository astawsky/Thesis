#!/usr/bin/env bash

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder


def put_all_graphs_into_a_big_grid(variables=phenotypic_variables):
    sns.set_context('paper')
    sns.set_style("ticks")
    
    # xticks = np.unique(output_df['distance'])[0::2] - 1
    # xticks = np.array(xticks, dtype=int)
    
    output_df['old_var'] = output_df['old_var'].replace(symbols['old_var'])
    output_df['young_var'] = output_df['young_var'].replace(symbols['young_var'])
    
    # top, bottom = output_df['correlations'].max(), output_df['correlations'].min()
    
    top, bottom = 1, -.5
    
    old_latex_symbols = {variable: symbols['old_var'][variable] for variable in variables}
    young_latex_symbols = {variable: symbols['young_var'][variable] for variable in variables}
    
    kinds = np.unique(output_df['dataset'])
    
    yticks = [-.5, -.25, 0, .25, .5, .75, 1]
    
    xmin, xmax = np.min(np.unique(output_df['distance'])), np.max(np.unique(output_df['distance']))

    xticks = np.array(np.unique(output_df['distance'])[0::2], dtype=int)
    
    fig, axes = plt.subplots(nrows=len(variables), ncols=len(variables), figsize=[7, 7])
    
    for row, row_var in zip(range(axes.shape[0]), list(old_latex_symbols.values())):
        for col, col_var in zip(range(axes.shape[1]), list(young_latex_symbols.values())):
            # define the axis we will be plotting on
            ax = axes[row, col]
            
            for kind in kinds:
                plot_df = output_df[(output_df['young_var'] == col_var) & (output_df['old_var'] == row_var) & (output_df['dataset'] == kind)]
                ax.plot(plot_df['distance'], plot_df['correlations'], marker='o', markersize=3, color=cmap[kind])
                
            # Setting the limits of all axes
            ax.set_ylim([bottom, top])
            ax.set_xlim([xmin, xmax])
            
            # Making the y-value = 0 easier to see
            ax.axhline(0, ls='-', color='black')
            
            # Adding the grids to the axes
            ax.yaxis.grid(True)
            ax.xaxis.grid(False)
            
            # To share the x,y-axes
            if row == axes.shape[0] - 1:
                ax.set_xlabel(col_var, size='x-large', loc='center')
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticks)
            else:
                ax.set_xlabel('')
                ax.set_xticklabels([])
            if col == 0:
                ax.set_yticks(yticks)
                ax.set_yticklabels(yticks)
                ax.set_ylabel(row_var+'        ', size='x-large', rotation=0, loc='center')
            else:
                ax.set_ylabel('')
                ax.set_yticks(yticks)
                ax.set_yticklabels([])
    
    plt.tight_layout(pad=.5)
    plt.savefig('{}/{}.png'.format(args.figs_location, args.inter_fig), dpi=300)
    # plt.show()
    plt.close()


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units_with_control.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages_with_control.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered_with_control.csv')
parser.add_argument('-cd', '--corr_dir', metavar='', type=str, help='Correlation dataframes in Data',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/total_inter_correlations.csv')
parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures')
parser.add_argument('-inter_fig', '--inter_fig', metavar='', type=str, help='What the big intergenerational figure name will be called.',
                    required=False, default='intergenerational_big_figure')

args = parser.parse_args()

create_folder(args.save_folder)

# Import the data
output_df = pd.read_csv(args.corr_dir)

# cut it at 10 generations
output_df = output_df[output_df['distance'] < 9]

cmap = {
    'physical units': 'blue',
    'trace centered': 'green'
}

# plot it
put_all_graphs_into_a_big_grid()
