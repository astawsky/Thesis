#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, cut_uneven_pairs, add_control, add_control_and_cut_extra_intervals, get_time_averages_df
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from itertools import combinations


def scatter_plot_three_datasets(tas):
    # Plot them for all combination of variables
    repeats = []
    for var1 in phenotypic_variables:
        for var2 in phenotypic_variables:
            # So we do not repeat the correlations
            if var2 in repeats:
                continue
            
            label = ''
            for dataset in tas.dataset.unique():
                corr = pearsonr(tas[(tas['dataset'] == dataset)].sort_values('trap_ID')[var1 + '_A'].values, tas[(tas['dataset'] == dataset)].sort_values('trap_ID')[var2 + '_B'].values)[0]
                
                label = label + '{}: {:.2},\n'.format(dataset, corr)
            
            label = label[:-2]
            
            seaborn_preamble()
            fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
            sns.scatterplot(data=tas, hue='dataset', x=var1 + '_A', y=var2 + '_B', hue_order=['SL', 'NL', 'CTRL'])
            ax.annotate(label, xy=(.5, .92), xycoords=ax.transAxes, fontsize=9, ha='center', va='bottom', color='red',
                        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
            plt.legend(title='')
            plt.xlabel(symbols['time_averages'][var1])
            plt.ylabel(symbols['time_averages'][var2])
            plt.show()
            plt.close()
        
        repeats.append(var1)
        
        
def scatter_plots_one_dataset(tas):
    # Create the folder where we will put the figures
    create_folder('{}/version 1'.format(args['data_origin']))
    
    # Plot them for all combination of variables
    repeats = []
    for var1 in phenotypic_variables:
        for var2 in phenotypic_variables:
            # So we do not repeat the correlations
            if var2 in repeats:
                continue
            
            # Discard the pairs that have an outlier
            relevant = tas.sort_values('trap_ID')[[var1 + '_A', var2 + '_B']].dropna()
            
            x, y = relevant[var1 + '_A'].values, relevant[var2 + '_B'].values
            
            # Calculate the two types of correlations
            corr = pearsonr(x, y)[0]
            spear = spearmanr(x, y)[0]
            
            seaborn_preamble()
            fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
            sns.scatterplot(data=tas, x=var1 + '_A', y=var2 + '_B')
            ax.annotate(r'$\rho=${:.2}, {:.2}, n={}'.format(corr, spear, len(relevant)), xy=(.5, .92), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
            plt.xlabel(symbols['time_averages'][var1] + r'$^A$')
            plt.ylabel(symbols['time_averages'][var2] + r'$^B$')
            plt.savefig('{}/version 1/{}_A and {}_B'.format(args['data_origin'], var1, var2), dpi=300)
            # plt.show()
            plt.close()
        
        repeats.append(var1)
        
        
def regression_plots(tas):
    # Create the folder where we will put the figures
    create_folder('{}/regression plots'.format(args['data_origin']))
    
    # Plot them for all combination of variables
    repeats = []
    for var1 in phenotypic_variables:
        for var2 in phenotypic_variables:
            # So we do not repeat the correlations
            if var2 in repeats:
                continue
            
            # Discard the pairs that have an outlier
            relevant = tas.sort_values('trap_ID')[[var1 + '_A', var2 + '_B']].dropna()
            
            x, y = relevant[var1 + '_A'].values, relevant[var2 + '_B'].values
            
            # Calculate the two types of correlations
            corr = pearsonr(x, y)[0]
            spear = spearmanr(x, y)[0]
            
            seaborn_preamble()
            fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
            sns.regplot(data=tas, x=var1 + '_A', y=var2 + '_B', line_kws={'color': 'black', 'ls': '--'})
            ax.annotate(r'$\rho=${:.2}, {:.2}, n={}'.format(corr, spear, len(relevant)), xy=(.5, .92), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
            plt.xlabel(symbols['time_averages'][var1] + r'$^A$')
            plt.ylabel(symbols['time_averages'][var2] + r'$^B$')
            plt.savefig('{}/regression plots/{}_A and {}_B'.format(args['data_origin'], var1, var2), dpi=300)
            # plt.show()
            plt.close()
        
        repeats.append(var1)


def main(args):
    pu = pd.read_csv(args['pu'])
    
    if 'CTRL' in args['input_args'].pair_lineages:
        pu = add_control(pu)
    
    # So we do not have any lineage-length effects between two pair lineages
    pu = cut_uneven_pairs(pu)
    
    # Where we will put the time-averages
    tas = pd.DataFrame(columns=[v + '_A' for v in phenotypic_variables] + [v + '_B' for v in phenotypic_variables] + ['trap_ID', 'dataset'])
    
    # For all the datasets separately
    for dataset in args['input_args'].pair_lineages:
        print(dataset)
        
        # For the trap IDs in this dataset
        for trap_id in pu[pu['dataset'] == dataset].trap_ID.unique():
            
            # The mean of the A/B traces
            a_to_add = pu[(pu['dataset'] == dataset) & (pu['trap_ID'] == trap_id) & (pu['trace'] == 'A')][phenotypic_variables].copy().mean().rename(
                dict(zip(phenotypic_variables, [v + '_A' for v in phenotypic_variables])))
            b_to_add = pu[(pu['dataset'] == dataset) & (pu['trap_ID'] == trap_id) & (pu['trace'] == 'B')][phenotypic_variables].copy().mean().rename(
                dict(zip(phenotypic_variables, [v + '_B' for v in phenotypic_variables])))
            
            # Put them together and set their Trap ID and dataset
            to_add = pd.concat([a_to_add, b_to_add])
            to_add['trap_ID'] = trap_id
            to_add['dataset'] = dataset
            
            # Add them to where we keep the time-averages/lineage-means
            tas = tas.append(to_add, ignore_index=True)
            
            # Make sure there are no NaNs
            if tas.isnull().values.any():
                print('error')
                print(tas[tas.isnull()])
                print()
                exit()
    
    # Sort them based on Trap ID
    tas = tas.sort_values('trap_ID')
    
    # # Save it
    # tas.to_csv('{}/tas.csv'.format(args['data_origin']), index=False)

    variables = [v + '_A' for v in phenotypic_variables] + [v + '_B' for v in phenotypic_variables]

    tas[variables] = tas[variables].where(
        np.abs(tas[variables] - tas[variables].mean()) < (3 * tas[variables].std()),
        other=np.nan
    )

    regression_plots(tas)
    scatter_plots_one_dataset(tas)
    
    var1, var2 = 'generationtime', 'growth_rate'
    
    # Discard the pairs that have an outlier
    relevant = tas.sort_values('trap_ID')[[var1 + '_A', var2 + '_B']].dropna()

    x, y = relevant[var1 + '_A'].values, relevant[var2 + '_B'].values

    # Calculate the two types of correlations
    corr = pearsonr(x, y)[0]
    spear = spearmanr(x, y)[0]
    
    fake_x = np.linspace(np.min(x), np.max(x))

    seaborn_preamble()
    fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
    # sns.regplot(data=tas, x=var1 + '_A', y=var2 + '_B', line_kws={'color': 'black', 'ls': '--'})
    sns.scatterplot(data=tas, x=var1 + '_A', y=var2 + '_B')
    plt.plot(fake_x, np.log(2) / fake_x, color='black', ls='--', label='log(2)')
    ax.annotate(r'$\rho=${:.2}, {:.2}, n={}'.format(corr, spear, len(relevant)), xy=(.5, .92), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    plt.xlabel(symbols['time_averages'][var1] + r'$^A$')
    plt.ylabel(symbols['time_averages'][var2] + r'$^B$')
    plt.savefig('{}/tau alpha ln2'.format(args['data_origin'], var1, var2), dpi=300)
    # plt.show()
    plt.close()


if __name__ == '__main__':
    import argparse
    import os
    
    # The variables we want to plot
    main_variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']
    
    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')
    
    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    parser.add_argument('-pair_lineages', '--pair_lineages', metavar='', nargs="+", help='Should the correlations be done for NL/SL/CTRL?', required=False,
                        default=['NL'])
    parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
                        default=dict(zip(['phenotypic_variables'], [phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    # Do all the Mother and Sister Machine data
    for data_origin in ['Pooled_SM']:  # ['Pooled_SM']:  # wang_datasets:  # input_args.dataset_names:
        print(data_origin)
        
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        create_folder(current_dir + '/' + data_origin)
        
        processed_data = os.path.dirname(current_dir) + '/Datasets/' + data_origin + '/ProcessedData/'
        
        """
        data_origin ==> Name of the dataset we are analysing
        raw_data ==> Where the folder containing the raw data for this dataset is
        processed_data ==> The folder we will put the processed data in
        """
        args = {
            'data_origin': data_origin,
            # Data singularities, long traces with significant filamentation, sudden drop-offs
            'save_figs': current_dir + '/' + data_origin,
            'pu': processed_data + 'z_score_under_3/physical_units_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'physical_units.csv',
            'input_args': input_args
        }
        
        main(args)
