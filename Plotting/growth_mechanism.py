#!/usr/bin/env bash

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, seaborn_preamble, mm_data_names


def main(args):
    
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
    
    # The variables we want to plot
    variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']
    
    # read the csv file where we keep the data
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    # read the csv file where we keep the data
    trace_centered, trace_means_df = pd.read_csv('{}/{}'.format(args.save_folder, args.tc)), pd.read_csv('{}/{}'.format(args.save_folder, args.ta))
    
    pyramid_of_pairwise_covariances()
    
    pyramid_of_pairwise_covariances(figurename='main variables without numbers', annot=False)
    
    variables = phenotypic_variables
    
    pyramid_of_pairwise_covariances(figsize=[7 * 1.5, 2.5 * 1.5], figurename='pyramids main variables')
    

if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    for data_origin in mm_data_names + ['SM']:
        print(data_origin)
        
        parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                            required=False, default='trace_centered.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        main(args)

        print('*'*200)
