#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, zscore
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, symbols, dataset_names
from string import ascii_uppercase as uppercase_letters


def main(args):
    def put_all_graphs_into_a_big_grid(df, label, variables=phenotypic_variables, remove_outliers=False, suffix='No suffix'):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
            
        print(df)

        for count, dataset in enumerate(df['dataset'].unique()):
            print(dataset)
            pair_a = df[(df['dataset'] == dataset) & (df['trace'] == 'A')].sort_values(['trap_ID', 'generation']).copy()[phenotypic_variables].rename(columns=latex_symbols)
            pair_b = df[(df['dataset'] == dataset) & (df['trace'] == 'B')].sort_values(['trap_ID', 'generation']).copy()[phenotypic_variables].rename(columns=latex_symbols)
    
            # Symmetrize them
            a_trace = pair_a.append(pair_b, ignore_index=True).reset_index(drop=True)
            b_trace = pair_b.append(pair_a, ignore_index=True).reset_index(drop=True)
            
            # take out the outliers
            if remove_outliers:
                before_a = len(pair_a)
                before_b = len(pair_b)
                pair_a = pair_a[(np.abs(zscore(pair_a, nan_policy='omit')) < 4).all(axis=1)].reset_index(drop=True)
                pair_b = pair_b[(np.abs(zscore(pair_b, nan_policy='omit')) < 4).all(axis=1)].reset_index(drop=True)
                after_a = len(pair_a)
                after_b = len(pair_b)
                print('percentage that stays after removing all values farther than the fourth standard deviation:', after_a / before_a, after_b / before_b)
            
            # now plot
            fig, axes = plt.subplots(nrows=len(variables) - 1, ncols=len(variables) - 1, figsize=[7, 7])
            
            for row, row_var in zip(range(axes.shape[0]), list(latex_symbols.values())[1:]):
                for col, col_var in zip(range(axes.shape[1]), list(latex_symbols.values())[:-1]):
    
                    # define the axis we will be plotting on
                    ax = axes[row, col]
                    
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
                    if col > row:
                        ax.set_frame_on(False)
                        ax.set_xticks([])
                        ax.set_yticks([])
                    else:
                        sns.regplot(x=a_trace[col_var], y=b_trace[row_var], ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
                        ax.annotate(str(np.corrcoef(a_trace[col_var], b_trace[row_var])[0, 1])[:4], xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
                    
                    if row == axes.shape[0] - 1:
                        ax.set_xlabel(col_var)
                    else:
                        ax.set_xlabel('')
                        ax.set_xticklabels([])
                    if col == 0:
                        ax.set_ylabel(row_var)
                    else:
                        ax.set_ylabel('')
                        ax.set_yticklabels([])
        
            plt.tight_layout(pad=.5)
            plt.savefig('{}/{}/{} {} {}.png'.format(args.figs_location, args.pc, label, dataset, suffix), dpi=300)
            # plt.show()
            plt.close()
    
    def heatmap_analogs(df, label, variables=phenotypic_variables, annot=False, suffix='No suffix'):
        
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}

        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[7, 2.5])
        axes = axes.flatten()
        
        for count, dataset in enumerate(df['dataset'].unique()):
            print(dataset)
            pair_a = df[(df['dataset'] == dataset) & (df['trace'] == 'A')].sort_values(['trap_ID', 'generation'])
            pair_b = df[(df['dataset'] == dataset) & (df['trace'] == 'B')].sort_values(['trap_ID', 'generation'])
            
            corr_df = pd.DataFrame(columns=variables, index=variables, dtype=float)
            
            # Symmetrize them
            a_trace = pair_a.append(pair_b, ignore_index=True).reset_index(drop=True)
            b_trace = pair_b.append(pair_a, ignore_index=True).reset_index(drop=True)
            
            repeats = []
            for var_a in variables:
                for var_b in variables:
                    if var_b not in repeats:
                        
                        corr_df[var_a].loc[var_b] = np.corrcoef(a_trace[var_a], b_trace[var_b])[0, 1]
                repeats.append(var_a)

            corr_df = corr_df[variables].rename(columns=latex_symbols, index=latex_symbols)

            mask = np.ones_like(corr_df)
            mask[np.tril_indices_from(mask)] = False

            if count == 2:
                cbar_ax = fig.add_axes([.91, .1, .03, .8])
                sns.heatmap(corr_df, annot=annot, square=True, vmin=-1, vmax=1, mask=mask, center=0, cbar=True, ax=axes[count], cbar_ax=cbar_ax, cbar_kws={"orientation": "vertical"}, fmt='.2f')
            else:
                sns.heatmap(corr_df, annot=annot, square=True, vmin=-1, vmax=1, mask=mask, center=0, cbar=False, ax=axes[count], fmt='.2f')
            axes[count].set_title(uppercase_letters[count], x=-.2, fontsize='xx-large')
            axes[count].set_xlabel(dataset)

        fig.tight_layout(rect=[0, 0, .9, 1])
        plt.savefig('{}/{}/{} {}.png'.format(args.figs_location, args.pch, label, suffix), dpi=300)
        plt.show()
        plt.close()
        
    # Create the folder where we will save the figures
    create_folder('{}/{}'.format(args.figs_location, args.pc))
    create_folder('{}/{}'.format(args.figs_location, args.pch))
    
    # import the labeled measured bacteria in trace-centered units
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.puc)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tcc)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    time_averages = pd.read_csv('{}/{}'.format(args.save_folder, args.tac)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    unique_ta = time_averages.drop_duplicates().sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])

    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        put_all_graphs_into_a_big_grid(df, kind, variables=phenotypic_variables, remove_outliers=True, suffix='full')

    # for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
    #     print(kind)
    #     put_all_graphs_into_a_big_grid(df, kind, variables=['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate'], remove_outliers=True, suffix='main five')

    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        heatmap_analogs(df, kind, variables=phenotypic_variables, annot=False, suffix='full')

    # for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
    #     print(kind)
    #     heatmap_analogs(df, kind, variables=['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate'], annot=True, suffix='main five')
    #
    # exit()

    # print('pu')
    # heatmap_analogs(physical_units, 'physical_units', variables=phenotypic_variables, annot=False)
    # print('tc')
    # heatmap_analogs(trace_centered, 'trace_centered', variables=phenotypic_variables, annot=False)
    # print('ta')
    # heatmap_analogs(time_averages, 'time_averages', variables=phenotypic_variables, annot=False)
    # print('ta unique')
    # heatmap_analogs(unique_ta, 'unique_ta', variables=phenotypic_variables, annot=False)
    #
    # print('pu')
    # put_all_graphs_into_a_big_grid(physical_units, 'physical_units', remove_outliers=True)
    # print('tc')
    # put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', remove_outliers=True)
    # print('ta')
    # put_all_graphs_into_a_big_grid(time_averages, 'time_averages', remove_outliers=True)
    # print('unique ta')
    # put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', remove_outliers=True)
    
    
if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    data_origin = 'SM'
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                        required=False, default='population_lineages.csv')
    parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                        required=False, default='time_averages.csv')
    parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                        required=False, default='ergodicity_breaking_parameter.csv')
    parser.add_argument('-kld', '--kld', metavar='', type=str,
                        help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                        required=False, default='kullback_leibler_divergences.csv')
    parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=False)
    args = parser.parse_args()
    
    main(args)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
