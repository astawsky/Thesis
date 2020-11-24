#!/usr/bin/env bash

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets, seaborn_preamble


def pyramid_of_pairwise_covariances(figsize=[7, 2.5], figurename='main variables without blanks', annot=True):
    
    def normalize_correlation(df, variables):
        cov = df.cov()
        corr_df = pd.DataFrame(columns=variables, index=variables, dtype=float)
        for param_1 in variables:
            for param_2 in variables:
                normalization = (physical_units[param_1].std() * physical_units[param_2].std())
                
                # # If the value lies in the noise range then make it zero
                # if -.1 < cov.loc[param_1, param_2] / normalization < .1:
                #     corr_df.loc[param_1, param_2] = float('nan')  # 0
                # else:
                #     corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
                corr_df.loc[param_1, param_2] = cov.loc[param_1, param_2] / normalization
        
        return corr_df
    
    seaborn_preamble()
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize)
    
    mask = np.ones_like(normalize_correlation(physical_units, variables))
    mask[np.tril_indices_from(mask)] = False
    vmax, vmin = 1, -1
    
    pu = normalize_correlation(physical_units, variables).rename(columns=symbols['physical_units'], index=symbols['physical_units'])
    ta = normalize_correlation(trace_means_df, variables).rename(columns=symbols['time_averages'], index=symbols['time_averages'])
    tc = normalize_correlation(trace_centered, variables).rename(columns=symbols['trace_centered'], index=symbols['trace_centered'])
    
    sns.heatmap(pu, annot=annot, center=0, vmax=vmax,
                vmin=vmin, cbar=False, ax=axes[0], mask=mask, square=True, fmt='.2f')
    axes[0].set_title('A', x=-.2, fontsize='xx-large')
    
    sns.heatmap(ta, annot=annot, center=0, vmax=vmax, vmin=vmin,
                cbar=False, ax=axes[1], mask=mask, square=True, fmt='.2f')
    axes[1].set_title('B', x=-.2, fontsize='xx-large')
    
    cbar_ax = fig.add_axes([.91, .1, .03, .8])
    
    sns.heatmap(tc, annot=annot, center=0, vmax=vmax,
                vmin=vmin, cbar=True, ax=axes[2], mask=mask, cbar_kws={"orientation": "vertical"}, square=True, cbar_ax=cbar_ax, fmt='.2f')
    axes[2].set_title('C', x=-.2, fontsize='xx-large')
    
    fig.tight_layout(rect=[0, 0, .9, 1])
    
    plt.savefig('{}/{}.png'.format(args.figs_location, figurename), dpi=300)
    # plt.show()
    plt.close()


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered.csv')
parser.add_argument('-bs', '--bootstraps', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)
parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures')


args = parser.parse_args()

create_folder(args.save_folder)

# The variables we want to plot
variables=['fold_growth', 'length_birth', 'generationtime', 'growth_rate']

# read the csv file where we keep the data
physical_units = pd.read_csv(args.pu)

# read the csv file where we keep the data
trace_centered, trace_means_df = pd.read_csv(args.tc), pd.read_csv(args.ta)

pyramid_of_pairwise_covariances()


pyramid_of_pairwise_covariances(figurename='main variables without numbers', annot=False)


variables = phenotypic_variables


pyramid_of_pairwise_covariances(figsize=[7*1.5, 2.5*1.5], figurename='pyramids all variables')
