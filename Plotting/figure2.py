#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.linear_model import LinearRegression
from string import ascii_uppercase as uppercase_letters
import argparse
            

""" Create the dataframe that contains the regression phenotypic_variables """


def slope_comparisons():
    # import the labeled measured bacteria in physical units
    cv_df = pd.read_csv('{}/variation_of_cumulative_time_averages.csv'.format(args.data_location))
    
    # This is where we will keep the regression phenotypic_variables
    slope_df = pd.DataFrame(columns=['Trace', 'Shuffled', 'Population', 'Trace-Centered'])
    intercept_df = slope_df.copy()
    
    # looping over the phenotypic_variables in highlights and their respective axes
    for count, param in enumerate(phenotypic_variables):
        to_add_intercept = {}
        to_add_slope = {}
        
        # Plot the four different kinds of lineages
        for label, marker in zip(slope_df.columns, ['o', 'x', '^', '+']):
            
            # the variance and generations of the TAs of this type of lineage
            dist = cv_df[(cv_df['param'] == param) & (cv_df['label'] == label)].sort_values('generation')
            
            # We get the best fit line
            lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'] + 1)).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
            
            # Add this regression phenotypic_variables to the dataframe
            to_add_slope.update({label: lin_reg.coef_[0]})
            to_add_intercept.update({label: lin_reg.intercept_[0]})
        
        # Add it to the dataframe to make a table or heatmap
        intercept_df = intercept_df.append(pd.DataFrame(to_add_intercept, index=[symbols['time_averages'][param]]))
        slope_df = slope_df.append(pd.DataFrame(to_add_slope, index=[symbols['time_averages'][param]]))
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style('ticks', {'axes.grid': True})

    fig, ax = plt.subplots(figsize=[7, 3.5])

    sns.heatmap(slope_df.T, annot=True, square=True, yticklabels=True)
    # plt.title('Slopes')
    ax.set_yticklabels(['Trace', 'Shuffled', 'Population', 'Trace-Centered'], rotation=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig('{}/Slopes heatmap.png'.format(args.figs_location), dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=[7, 3.5])

    sns.heatmap(intercept_df.T, annot=True, square=True, yticklabels=True)
    # plt.title('Intercepts')
    ax.set_yticklabels(['Trace', 'Shuffled', 'Population', 'Trace-Centered'], rotation=0)
    plt.tight_layout()
    # plt.show()
    plt.savefig('{}/Intercepts heatmap.png'.format(args.figs_location), dpi=300)
    plt.close()
            
            
""" This creates the big figure 2 from the paper. Includes the individual figures code where  """


def variance_timeaverage_per_generation():
    # import the labeled measured bacteria in physical units
    cv_df = pd.read_csv('{}/variation_of_cumulative_time_averages.csv'.format(args.data_location))
    
    # The cycle phenotypic_variables we will look at
    highlights = ['length_birth', 'growth_rate', 'fold_growth']
    
    # set a style on seaborn module
    sns.set_context('paper')
    sns.set_style('ticks', {'axes.grid': True})
    
    # create the figure and construct the layout of the figure
    fig, axes = plt.subplots(nrows=1, ncols=len(highlights), sharey=True, figsize=[7, 3])
    
    # looping over the phenotypic_variables in highlights and their respective axes
    for count, param in enumerate(highlights):
        # Specify the axis
        ax = axes[count]
        
        # Plot the four different kinds of lineages
        for label, marker in zip(['Trace', 'Population', 'Trace-Centered', 'Shuffled'], ['o', 'x', '^', '+']):
            
            # the variance and generations of the TAs of this type of lineage
            dist = cv_df[(cv_df['param'] == param) & (cv_df['label'] == label)].sort_values('generation')
            
            # We get the best fit line
            lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'] + 1)).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
            
            # plot the plots
            ax.scatter(dist['generation'] + 1, dist['var'], marker=marker, alpha=.5)
            ax.plot(dist['generation'] + 1, np.exp(lin_reg.predict(np.array(np.log(dist['generation'] + 1)).reshape(-1, 1)).flatten()), ls='--', label='')
            
            # reg_order = PowerRegression(dist['generation'] + 1, dist['var'], bs=500)
            # ax.scatter(dist['generation'] + 1, dist['var'], label=reg_order.param_labels['a'], marker=marker, alpha=.5)
            # ax.plot(dist['generation'] + 1, reg_order.prediction, ls='--', label='')

        # Add the grid
        ax.grid(True)
        # ax.legend()
        ax.set(xscale="log", yscale="log")
        # ax.tick_params(direction='out')
        ax.set_xlabel(r'$N$')
        ax.set_ylabel(r'$Var(${}$)$'.format(symbols['time_averages'][param]))
        # ax.set_ylim([low_bound, high_bound])
        # ax.set_xlim(left=0)

        ax.set_title(uppercase_letters[count], x=-.15, fontsize='xx-large')
        
    # Save the big figure 2.png
    fig.tight_layout()
    # plt.show()
    plt.savefig('{}/Figure 2.png'.format(args.figs_location), dpi=300)
    plt.close()
    

parser = argparse.ArgumentParser(description='Plot the dependance of the time-averages on the lineage length.')
parser.add_argument('-d', '--data_location', metavar='', type=str, help='Where the dataframes are saved.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures')


args = parser.parse_args()


create_folder(args.data_location)
    
# phenotypic_variables, symbols, args
variance_timeaverage_per_generation()


slope_comparisons()
