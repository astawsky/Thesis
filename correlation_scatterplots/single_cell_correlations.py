#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import linregress, pearsonr
import pingouin as pg
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from AnalysisCode.global_variables import (
    phenotypic_variables, create_folder, symbols, units, dataset_names, get_time_averages_df, sm_datasets, wang_datasets, tanouchi_datasets, seaborn_preamble,
    cgsc_6300_wang_exps, cmap, pool_experiments, lexA3_wang_exps, mm_datasets
)


def main_groups(pu, ta, tc, name):
    def put_all_graphs_into_a_big_grid(df, label, filename, variables=phenotypic_variables):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        unit_symbols = {variable: units[variable] if label != 'trace_centered' else '' for variable in variables}
        
        # Take out the outliers
        df[variables] = df[variables].where(
            np.abs(df[variables] - df[variables].mean()) < (3 * df[variables].std()),
            other=np.nan
        )
        
        # Replace the phenotypic variables with their latex counterparts
        df = df[variables].copy().rename(columns=latex_symbols)
        
        # In case we just want to see one scatterplot
        fig, axes = plt.subplots(nrows=len(variables) - 1, ncols=len(variables) - 1, figsize=[7, 7])
        
        for row, row_var in zip(range(axes.shape[0]), list(latex_symbols.values())[1:]):
                for col, col_var in zip(range(axes.shape[1]), list(latex_symbols.values())[:-1]):
                    # define the axis we will be plotting on
                    ax = axes[row, col]
                    
                    # So it looks presentable
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
                    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    
                    # Lower diagonal matrix
                    if col > row:
                        ax.set_frame_on(False)
                        ax.set_xticks([])
                        ax.set_yticks([])
                    else:
                        index_without_nans = [index for index in df.index if ~np.isnan(df[[col_var, row_var]].loc[index].values).any()]
                        
                        col_array = df[col_var].loc[index_without_nans].values
                        row_array = df[row_var].loc[index_without_nans].values
                        
                        pcorr = str(pearsonr(col_array, row_array)[0])[:4]
                        slope, intercept, _, _, std_err = linregress(col_array, row_array)
                        
                        # pcorr = str(pg.corr(df[col_var].values, df[row_var].values, method='pearson')['r'].loc['pearson'].round(2))[:4]
                        # slope, intercept = str(linregress(df[list(latex_symbols.values())[1]], df[list(latex_symbols.values())[0]])[0])[:4]
                        
                        # sns.regplot(x=df[col_var], y=df[row_var], data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'color': 'grey'})
                        # ax.scatter(df[col_var].values, df[row_var].values, color='gray')
                        # ax.plot(df[col_var].values, [intercept + slope * vel for vel in df[col_var].values])
                        
                        ax.scatter(col_array, row_array, color=cmap[0], marker='o', alpha=.3)  # .1
                        ax.plot(np.unique(col_array), [intercept + slope * vel for vel in np.unique(col_array)], c='k', ls='--')
                        ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {:.2}$'.format(np.round(slope, 2)), xy=(.5, .8), xycoords=ax.transAxes, fontsize=7, ha='center', va='bottom', color='red',
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
        # In case we just want to see one scatterplot
        plt.savefig(filename, dpi=300)
        plt.show()
        plt.close()

    print('unique ta')
    put_all_graphs_into_a_big_grid(ta, 'time_averages', f'/Users/alestawsky/PycharmProjects/Thesis/correlation_scatterplots/{name} ta.png',
                                   variables=phenotypic_variables)


def main(args):
    def put_all_graphs_into_a_big_grid(df, label, variables=phenotypic_variables, suffix='', reg_dict=None):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        unit_symbols = {variable: units[variable] if label != 'trace_centered' else '' for variable in variables}
        
        # Take out the outliers
        df[variables] = df[variables].where(
            np.abs(df[variables] - df[variables].mean()) < (3 * df[variables].std()),
            other=np.nan
        )
        
        # Replace the phenotypic variables with their latex counterparts
        df = df[variables].copy().rename(columns=latex_symbols)
        
        # In case we just want to see one scatterplot
        if len(variables) == 2:
            
            sns.set_context('talk')
            
            fig, axes = plt.subplots(figsize=[5, 5], tight_layout=True)
            
            sym1 = list(latex_symbols.values())[0]
            sym2 = list(latex_symbols.values())[1]
            
            relevant = df[[sym2, sym1]].dropna()
            
            # if len(np.unique(variables)) == 1:
            #     sym1 = list(latex_symbols.values())[0]
            #     sym2 = list(latex_symbols.values())[0]
            # else:
            #     sym1 = list(latex_symbols.values())[0]
            #     sym2 = list(latex_symbols.values())[1]
            #
            # print(list(latex_symbols.values()))
            # print(df[sym1].values)
            # print(df[sym2].values)
            
            pcorr = str(pearsonr(relevant[sym2].values, relevant[sym1].values)[0])[:4]
            # pcorr = str(pg.corr(df[sym2], df[sym1], method='pearson')['r'].loc['pearson'].round(2))[:4]
            slope, intercept, r_value, _, std_err = linregress(relevant[sym2], relevant[sym1])
            
            reg_dict.update({
                len(reg_dict): {
                    'kind': label,
                    'slope': slope,
                    'intercept': intercept,
                    'std_err': std_err,
                    'r_value': r_value
                }
            })
            
            # axes.scatter(relevant[sym2], relevant[sym1], color='gray', marker='o', alpha=.1)
            # # sns.kdeplot(x=relevant[sym2], y=relevant[sym1], ax=axes, color='grey')
            # axes.plot(np.unique(relevant[sym2]), [intercept + slope * vel for vel in np.unique(relevant[sym2])])
            
            sns.regplot(data=relevant, x=sym2, y=sym1, line_kws={'color': 'black', 'ls': '--', 'lw': 2})
            
            # sns.regplot(x=df[sym2], y=df[sym1], data=df, ax=axes, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
            
            axes.annotate(r'$\rho = $' + pcorr + ', ' + r'$\beta = {}$'.format(str(slope)[:4]), xy=(.5, .72), xycoords=axes.transAxes, fontsize=13, ha='center', va='bottom',
                          color='red',
                          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
            axes.set_ylabel(sym1 + ' ' + list(unit_symbols.values())[0])
            axes.set_xlabel(sym2 + ' ' + list(unit_symbols.values())[1])
            # # plt.tight_layout(pad=.5)
            # plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scc, label, suffix), dpi=300)
            # # plt.show()
            # plt.close()
        else:  # If we want the matrix of scatterplots
            
            fig, axes = plt.subplots(nrows=len(variables) - 1, ncols=len(variables) - 1, figsize=[7, 7])
            
            for row, row_var in zip(range(axes.shape[0]), list(latex_symbols.values())[1:]):
                for col, col_var in zip(range(axes.shape[1]), list(latex_symbols.values())[:-1]):
                    # define the axis we will be plotting on
                    ax = axes[row, col]
                    
                    # So it looks presentable
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2g'))
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
                    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                    
                    # Lower diagonal matrix
                    if col > row:
                        ax.set_frame_on(False)
                        ax.set_xticks([])
                        ax.set_yticks([])
                    else:
                        index_without_nans = [index for index in df.index if ~np.isnan(df[[col_var, row_var]].loc[index].values).any()]
                        
                        col_array = df[col_var].loc[index_without_nans].values
                        row_array = df[row_var].loc[index_without_nans].values
                        
                        pcorr = str(pearsonr(col_array, row_array)[0])[:4]
                        slope, intercept, _, _, std_err = linregress(col_array, row_array)
                        
                        # pcorr = str(pg.corr(df[col_var].values, df[row_var].values, method='pearson')['r'].loc['pearson'].round(2))[:4]
                        # slope, intercept = str(linregress(df[list(latex_symbols.values())[1]], df[list(latex_symbols.values())[0]])[0])[:4]
                        
                        # sns.regplot(x=df[col_var], y=df[row_var], data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'color': 'grey'})
                        # ax.scatter(df[col_var].values, df[row_var].values, color='gray')
                        # ax.plot(df[col_var].values, [intercept + slope * vel for vel in df[col_var].values])
                        
                        ax.scatter(col_array, row_array, color=cmap[0], marker='o', alpha=.3)  # .1
                        ax.plot(np.unique(col_array), [intercept + slope * vel for vel in np.unique(col_array)], c='k', ls='--')
                        ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {:.2}$'.format(np.round(slope, 2)), xy=(.5, .8), xycoords=ax.transAxes, fontsize=7, ha='center', va='bottom', color='red',
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
        # In case we just want to see one scatterplot
        if len(variables) == 2:
            # plt.savefig(args['Figures'] + '/{}/separate/{}/'.format(label, args['data_origin']) + '{}_{}{}.png'.format(variables[0], variables[1], suffix), dpi=300)
            plt.show()
            plt.close()
            
            return reg_dict
        else:
            # plt.savefig(args['Figures'] + '/{}/together/'.format(label) + '{}{}.png'.format(args['data_origin'], suffix), dpi=300)
            plt.show()
            plt.close()
    
    # for the regression parameters
    regression = {}
    
    # import the labeled measured bacteria in trace-centered units
    physical_units = pd.read_csv(args['pu'])
    trace_centered = pd.read_csv(args['tc'])
    
    # Calculate this
    time_averages = get_time_averages_df(physical_units, phenotypic_variables)[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
    time_averages = time_averages[time_averages['max_gen'] > 7]
    
    # for kind, dataframe in zip(['physical_units', 'trace_centered', 'time_averages'], [physical_units, trace_centered, time_averages]):
    #     print(kind)
    #
    #     # # Plot all the variables -- Heatmap
    #     # heatmap_analogs(dataframe, kind, variables=phenotypic_variables, suffix='')
    #
    #     # # Plot all the variables -- Scatter Regression Plot
    #     # put_all_graphs_into_a_big_grid(dataframe, kind, variables=phenotypic_variables, suffix='')
    #
    #     # Scatter Regression Plots for each variable pair individually
    #     repeats = []
    #     for var1 in phenotypic_variables:
    #         for var2 in phenotypic_variables:
    #             if var2 in repeats or var2 == var1:
    #                 continue
    #
    #             var1, var2 = 'generationtime', 'growth_rate'
    #
    #             # the separate figures
    #             regression = put_all_graphs_into_a_big_grid(dataframe, kind, variables=[var1, var2], suffix='', reg_dict=regression)
    #
    #         # So that we don't repeat the scatter plots
    #         repeats.append(var1)
    #
    # # regression dataframe for analysis
    # pd.DataFrame.from_dict(regression, "index").to_csv('regression_dataframes/'+args['data_origin']+'_regressions.csv', index=False)
    
    # # Plot all the variables -- Heatmap
    # print('pu')
    # heatmap_analogs(physical_units, 'physical_units', variables=phenotypic_variables, suffix='')
    # print('tc')
    # heatmap_analogs(trace_centered, 'trace_centered', variables=phenotypic_variables, suffix='')
    # print('ta unique')
    # heatmap_analogs(unique_ta, 'time_averages', variables=phenotypic_variables, suffix='')
    #
    # Plot all the variables -- Scatter Regression Plot
    # print('pu')
    # put_all_graphs_into_a_big_grid(physical_units, 'physical_units', variables=phenotypic_variables, suffix='')
    # print('tc')
    # put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', variables=phenotypic_variables, suffix='')
    print('unique ta')
    put_all_graphs_into_a_big_grid(time_averages, 'time_averages', variables=phenotypic_variables, suffix='')
    #
    # # Scatter Regression Plots for each variable pair individually
    # repeats = []
    # for var1 in phenotypic_variables[3:]:
    #     for var2 in phenotypic_variables[3:]:
    #         if var2 in repeats:
    #             continue
    #         print(var1, var2)
    #
    #         # the separate figures
    #         put_all_graphs_into_a_big_grid(physical_units, 'physical_units', variables=[var1, var2], suffix='')
    #     repeats.append(var1)


if __name__ == '__main__':
    import argparse
    import os
    
    # The variables we want to plot
    main_variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']
    
    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')
    
    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    
    # Finalize the arguments
    input_args = parser.parse_args()

    for name, group in zip(['Wang 2010, CGSC 6300', 'Wang 2010, lexA3', 'Susman and Kohram', 'Vashistha 2021', 'Tanouchi 25C', 'Tanouchi 27C', 'Tanouchi 37C'],
                           [cgsc_6300_wang_exps, lexA3_wang_exps, mm_datasets, ['Pooled_SM'], ['MC4100_25C (Tanouchi 2015)'], ['MC4100_27C (Tanouchi 2015)'], ['MC4100_37C (Tanouchi 2015)']]):
        print(name, group)
        
        pu = pool_experiments(group, name, dimensions='pu')
        ta = pool_experiments(group, name, dimensions='ta')[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
        ta = ta[ta['max_gen'] > 7]
        tc = pool_experiments(group, name, dimensions='tc')
        
        main_groups(pu, ta, tc, name)
        
    
    # Do all the Mother and Sister Machine data
    for data_origin in ['Pooled_SM'] + mm_data + tanouchi_datasets + cgsc_6300_wang_exps:  # input_args.dataset_names:
        print(data_origin)
        
        filepath = os.path.dirname(os.path.abspath(__file__))
        
        create_folder(filepath + '/regression_dataframes')
        
        processed_data = os.path.dirname(filepath) + '/Datasets/' + data_origin + '/ProcessedData/'
        
        """
        data_origin ==> Name of the dataset we are analysing
        raw_data ==> Where the folder containing the raw data for this dataset is
        processed_data ==> The folder we will put the processed data in
        """
        args = {
            'data_origin': data_origin,
            'MM': False if data_origin in sm_datasets else True,
            # Data singularities, long traces with significant filamentation, sudden drop-offs
            'Figures': filepath + '/Figures',
            'pu': processed_data + 'z_score_under_3/physical_units_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'physical_units.csv',
            'tc': processed_data + 'z_score_under_3/trace_centered_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'trace_centered.csv'
        }
        
        # Where we will put the figures
        create_folder(args['Figures'])
        for lll in ['physical_units', 'trace_centered', 'time_averages']:
            create_folder(args['Figures'] + '/{}'.format(lll))
            create_folder(args['Figures'] + '/{}/together'.format(lll))
            create_folder(args['Figures'] + '/{}/separate'.format(lll))
            create_folder(args['Figures'] + '/{}/separate/{}'.format(lll, data_origin))
        
        main(args)
        
        print('*' * 200)
