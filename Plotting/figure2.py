#!/usr/bin/env bash

import sys

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, seaborn_preamble
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.linear_model import LinearRegression
from string import ascii_uppercase as uppercase_letters
import argparse
            

""" Create the dataframe that contains the regression phenotypic_variables """


def slope_comparisons(args, variation_of_expanding_mean):
    
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
            dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
            
            # So each kind of lineage starts at the 1 and keeps their slope
            dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
            
            # We get the best fit line
            lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
            
            # Add this regression phenotypic_variables to the dataframe
            to_add_slope.update({label: lin_reg.coef_[0]})
            to_add_intercept.update({label: lin_reg.intercept_[0]})
        
        # Add it to the dataframe to make a table or heatmap
        intercept_df = intercept_df.append(pd.DataFrame(to_add_intercept, index=[symbols['time_averages'][param]]))
        slope_df = slope_df.append(pd.DataFrame(to_add_slope, index=[symbols['time_averages'][param]]))
    
    # set a style on seaborn module
    seaborn_preamble()

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


def variance_timeaverage_per_generation(args, variation_of_expanding_mean):
    
    # The cycle phenotypic_variables we will look at
    highlights = ['length_birth', 'growth_rate', 'generationtime']
    
    # set a style on seaborn module
    seaborn_preamble()
    
    # create the figure and construct the layout of the figure
    fig, axes = plt.subplots(nrows=1, ncols=len(highlights), sharey=True, figsize=[7, 3])
    
    # looping over the phenotypic_variables in highlights and their respective axes
    for count, param in enumerate(highlights):
        # Specify the axis
        ax = axes[count]
        
        # Plot the four different kinds of lineages
        for label, marker in zip(['Trace', 'Population', 'Trace-Centered', 'Shuffled'], ['o', 'x', '^', '+']):
            
            # the variance and generations of the TAs of this type of lineage
            dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
            
            # So each kind of lineage starts at the 1 and keeps their slope
            dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
            
            # # We get the best fit line
            # lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
            # ax.plot(dist['generation'], np.exp(lin_reg.predict(np.array(np.log(dist['generation'])).reshape(-1, 1)).flatten()), ls='--', label='')
            
            # plot the plots
            ax.scatter(dist['generation'], dist['var'], marker=marker, alpha=.5)
            
            # reg_order = PowerRegression(dist['generation'], dist['var'], bs=500)
            # ax.scatter(dist['generation'], dist['var'], label=reg_order.param_labels['a'], marker=marker, alpha=.5)
            # ax.plot(dist['generation'], reg_order.prediction, ls='--', label='')

        # Add the grid
        ax.grid(True)
        # ax.legend()
        ax.set(xscale="log", yscale="log")
        # ax.tick_params(direction='out')
        ax.set_xlabel(r'$N$')
        ax.set_ylabel(r'$Var(${}$)$'.format(symbols['time_averages'][param], symbols['physical_units'][param]))
        # ax.set_ylim([low_bound, high_bound])
        # ax.set_xlim(left=0)

        ax.set_title(uppercase_letters[count], x=-.15, fontsize='xx-large')
        
    # Save the big figure 2.png
    fig.tight_layout()
    # plt.show()
    plt.savefig('{}/Figure 2.png'.format(args.figs_location), dpi=300)
    plt.close()
    
    
""" This creates the rest of the figures left out of fig 2 in the main text to be in the appendix """


def variance_timeaverage_per_generation_rest_of_figures(args, variation_of_expanding_mean):
    
    # The cycle phenotypic_variables we will look at
    highlights = ['fold_growth', 'length_final', 'division_ratio', 'added_length']
    
    # set a style on seaborn module
    seaborn_preamble()
    
    # create the figure and construct the layout of the figure
    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, figsize=[7, 7])
    axes = axes.flatten()
    
    # looping over the phenotypic_variables in highlights and their respective axes
    for count, param in enumerate(highlights):
        # Specify the axis
        ax = axes[count]
        
        # Plot the four different kinds of lineages
        for label, marker in zip(['Trace', 'Population', 'Trace-Centered', 'Shuffled'], ['o', 'x', '^', '+']):
            
            # the variance and generations of the TAs of this type of lineage
            dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == param) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')
            
            # So each kind of lineage starts at the 1 and keeps their slope
            dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values

            # # We get the best fit line
            # lin_reg = LinearRegression(fit_intercept=True).fit(np.array(np.log(dist['generation'])).reshape(-1, 1), np.array(np.log(dist['var'])).reshape(-1, 1))
            # ax.plot(dist['generation'], np.exp(lin_reg.predict(np.array(np.log(dist['generation'])).reshape(-1, 1)).flatten()), ls='--', label='')
            
            # plot the plots
            ax.scatter(dist['generation'], dist['var'], marker=marker, alpha=.5)
            
            # reg_order = PowerRegression(dist['generation'], dist['var'], bs=500)
            # ax.scatter(dist['generation'], dist['var'], label=reg_order.param_labels['a'], marker=marker, alpha=.5)
            # ax.plot(dist['generation'], reg_order.prediction, ls='--', label='')

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
    plt.savefig('{}/Rest of Figure 2.png'.format(args.figs_location), dpi=300)
    plt.close()
    

def main(args):
    # import the labeled measured bacteria in physical units
    variation_of_expanding_mean = pd.read_csv('{}/variation_of_expanding_mean.csv'.format(args.save_folder))
    
    variance_timeaverage_per_generation_rest_of_figures(args, variation_of_expanding_mean)
    
    variance_timeaverage_per_generation(args, variation_of_expanding_mean)
    
    slope_comparisons(args, variation_of_expanding_mean)
