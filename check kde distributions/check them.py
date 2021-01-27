#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations


def check(chosen_datasets, save_fig, lin_type):

    seaborn_preamble()
    fig, ax = plt.subplots(tight_layout=True, figsize=[5, 5])
    total_df = pd.DataFrame()
    count = 0
    for data_origin in chosen_datasets:
        print(data_origin)
        pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(pu)
        pu['experiment'] = data_origin
        
        total_df = total_df.append(pu, ignore_index=True)
    
        sns.kdeplot(data=pu, x='growth_rate', y='generationtime', label=data_origin, color=cmap[count], alpha=.5)
        count += 1
        plt.scatter(pu['growth_rate'].mean(), pu['generationtime'].mean(), color=cmap[count], alpha=.8)
    
    plt.plot(np.linspace(.5, 2), np.log(2) / np.linspace(.5, 2), ls='--', c='k', alpha=.7)
    # plt.ylim([0, 1.4])
    plt.xlim([.25, 2.8])
    plt.xlabel(symbols['physical_units']['growth_rate'])
    plt.ylabel(symbols['physical_units']['generationtime'])
    plt.legend()
    plt.show()
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
    parser.add_argument('-kinds_of_correlations', '--kinds_of_correlations', metavar='', nargs="+", help='Calculate pearson and/or variance decomposition correlation?', required=False,
                        default=['decomposition', 'pearson'])
    parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
                        default=dict(zip(['phenotypic_variables'], [phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    check(sm_datasets[:-1], os.path.dirname(os.path.abspath(__file__)), lin_type='datasets')  # mm_datasets +
