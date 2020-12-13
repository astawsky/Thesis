#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import linregress, zscore
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, symbols, limit_lineage_length, trace_center_a_dataframe
from string import ascii_uppercase as uppercase_letters


def main():
    
    def shuffle_and_recreate(trajectory):
        
        new_lineages = {'Trace': trajectory['length_birth'].values}
        
        for what_to_shuffle in ['fold_growth', 'divison_ratio', 'fold_growth and division_ratio']:
            x = [trajectory[trajectory['generation'] == 0]['length_birth']]
            print(what_to_shuffle)
            if what_to_shuffle == 'fold_growth':
                fg = trajectory[what_to_shuffle].sample(frac=1, replace=False).values
                dr = trajectory['division_ratio'].values
            elif what_to_shuffle == 'division_ratio':
                fg = trajectory['fold_growth'].values
                dr = trajectory[what_to_shuffle].sample(frac=1, replace=False).values
            else:
                fg = trajectory['fold_growth'].sample(frac=1, replace=False).values
                dr = trajectory['division_ratio'].sample(frac=1, replace=False).values
                
            for d, f in zip(dr[:-1], fg[:-1]):
                x.append(x[-1]*d*np.exp(f))
            
            new_lineages.update({'shuffled '+what_to_shuffle: x})
            
        return new_lineages
    
    variable='length_birth'
    
    # import the data
    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/physical_units.csv')

    apply_mask = lambda df: df[(df['trap_ID'] == trap_id) & (df['dataset'] == dataset) & (df['trace'] == trace)].sort_values('generation').copy()

    # plot each lineage
    for dataset in physical_units['dataset'].unique():
    
        for trap_id in physical_units[(physical_units['dataset'] == dataset)]['trap_ID'].unique():
        
            for trace in ['A', 'B']:
            
                lineage = apply_mask(physical_units)
                
                if len(lineage) > 20:
                    
                    new_lineages = shuffle_and_recreate(lineage)
                    
                    for k,v in new_lineages.items():
                        print(k, v)
                        print(len(v))
                        plt.plot(v, marker='.', label=k)
                    plt.legend()
                    plt.show()
                    plt.close()
    
    
    
    
    
    
    
    
    
    
    # physical_units = limit_lineage_length(physical_units, min_gens=50)
    # trace_centered = trace_center_a_dataframe(physical_units)
    
    # trace_centered = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/trace_centered.csv')
    population_sampled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/population_lineages.csv')
    # shuffled_pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/shuffled_generations.csv')
    # tc_shuffled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/tc_shuffled.csv')
    variation_of_expanding_mean = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/variation_of_expanding_mean.csv')
    
    kinds = ['Trace', 'Population']
    dataframes = [physical_units, population_sampled]
    
    # fig, axes = plt.subplots(nrows=4, ncols=3, tight_layout=True, sharex='col')
    fig, axes = plt.subplots(nrows=4, ncols=3, tight_layout=True, figsize=[3 * 2, 4 * 2], sharex='col')
    
    for ax in axes[:, 1]:
        ax.set_yticklabels([])
    
    for col, kind in zip(np.arange(axes.shape[1]), kinds):
        axes[0, col].set_title(kind)
        if col != axes.shape[1] - 1:
            axes[-1, col].set_xlabel(r'$N$')
        else:
            axes[-1, col].set_xlabel(r'$log(N)$')
    
    for row, variable in enumerate(['fold_growth', 'generationtime', 'length_birth', 'growth_rate']):
        axes[row, 0].set_ylabel(symbols['time_averages'][variable] + r'$_{N}$' if variable not in ['length_birth', 'length_final'] else symbols['time_averages'][variable] + r'$_{,N}$')
        
        for col, kind in zip(np.arange(axes.shape[1] - 1), kinds):
            axes[row, col].set_ylim([min([df[df['generation'] == 0][variable].min() for df in dataframes]), max([df[df['generation'] == 0][variable].max() for df in dataframes])])
        
        # plot each lineage
        for dataset in physical_units['dataset'].unique():
            
            for trap_id in physical_units[(physical_units['dataset'] == dataset)]['trap_ID'].unique():
                
                for trace in ['A', 'B']:
                    
                    apply_mask = lambda df: df[(df['trap_ID'] == trap_id) & (df['dataset'] == dataset) & (df['trace'] == trace)][variable].copy().expanding().mean().values
                    
                    for df, ax in zip(dataframes, axes[row, :-1]):
                        lineage = apply_mask(df)
                        
                        # start plotting
                        ax.plot(np.arange(1, 1 + len(lineage)), lineage, marker=',')
        
        # show the variance analysis for the 5 different kinds of lineages
        for kind, marker, color in zip(kinds, ['.', 'x'], sns.color_palette('tab10')):
            
            # the variance and generations of the TAs of this type of lineage
            dist = variation_of_expanding_mean[
                (variation_of_expanding_mean['param'] == variable) & (variation_of_expanding_mean['label'] == kind) & (variation_of_expanding_mean['generation'] <= 49)].sort_values('generation')
            
            # print(dist)
            
            # So each kind of lineage starts at the 1 and keeps their slope
            # dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
            
            # print(dist)
            # exit()
            
            axes[row, -1].scatter(dist['generation'], dist['var'], color=color, marker=marker, alpha=.7)
            if kind == 'Population':
                slope, intercept, _, _, std_err = linregress(np.log(dist['generation'].values), np.log(dist['var'].values))
                axes[row, -1].plot(dist['generation'], np.exp([intercept + slope * np.log(gen) for gen in dist['generation']]), ls='--', color=color)
                axes[row, -1].plot(dist['generation'], np.exp([intercept - 1 * np.log(gen) for gen in dist['generation']]), ls='-.', color='black')
            
            # Add the grid
            axes[row, -1].grid(True)
            # ax.legend()
            axes[row, -1].set(xscale="log", yscale="log")
            # ax.tick_params(direction='out')
            axes[row, -1].set_xlabel(r'$log(N)$')
            axes[row, -1].set_ylabel(
                r'$Var($' + symbols['time_averages'][variable] + r'$_{N})$' if variable not in ['length_birth', 'length_final'] else r'$Var($' + symbols['time_averages'][variable] + r'$_{,N})$')
            # axes[row, -1].set_ylabel(
            #     r'$Var($' + symbols['time_averages'][variable] + r'$_{N}) / Var($' + symbols['time_averages'][variable] + '$_{0})$' if variable not in ['length_birth', 'length_final'] else
            #     r'$Var($' + symbols['time_averages'][variable] + r'$_{,N}) / Var($' + symbols['time_averages'][variable] + '$_{,0})$')
        
        # Set the titles saying what kind of lineage we are looking at
        # and also show the mean of all the lineages across all generation-time
        for ax, df, kind in zip(axes[row, :-1], dataframes, kinds):
            ax.axhline(df[variable].mean(), ls='-', color='black')
            ax.set_xlim([1, 50])
            ax.yaxis.grid(True)
            # ax.set_xlabel(r'$N$')
            # ax.set_ylabel(symbols['time_averages'][variable] + r'$_{N}$' if variable not in ['length_birth', 'length_final'] else symbols['time_averages'][variable] + r'$_{,N}$')
    
    plt.savefig('Better figure 2.png', dpi=300)
    # plt.show()
    plt.close()


if __name__ == '__main__':
    main()
