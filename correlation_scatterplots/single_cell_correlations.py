#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import linregress, pearsonr
import pingouin as pg
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from AnalysisCode.global_variables import phenotypic_variables, create_folder, symbols, units, dataset_names, get_time_averages_df, sm_datasets, wang_datasets, tanouchi_datasets


def previous_division_scatterplots(df, label, var_prev, var_after, suffix=''):
    prev_array = np.concatenate([
        df[(df['lineage_ID'] == lineage_id)].sort_values('generation')[var_prev].values[:-1] for lineage_id in df.lineage_ID.unique()]).flatten()
    after_array = np.concatenate([
        df[(df['lineage_ID'] == lineage_id)].sort_values('generation')[var_after].values[1:] for lineage_id in df.lineage_ID.unique()]).flatten()
    
    # prev_array = prev_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]
    # after_array = after_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]
    
    # # Remove outliers
    # for ar in [prev_array, after_array]:
    #     before = len(ar)
    #     ar = ar[np.abs(zscore(ar)) < 3]
    #     after = len(ar)
    #     print('percentage that stays after removing all values farther than the fourth standard deviation:', after / before)
    
    pcorr = str(pg.corr(prev_array, after_array, method='spearman')['r'].loc['spearman'])[:4]
    slope = str(linregress(prev_array, after_array)[0])[:4]
    
    prev_units = units[var_prev] if label != 'trace_centered' else ''
    after_units = units[var_after] if label != 'trace_centered' else ''
    
    sns.set_context('talk')
    
    fig, ax = plt.subplots(figsize=[3.5, 3.5], tight_layout=True)
    
    sns.regplot(x=prev_array, y=after_array, data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
    ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(slope), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    
    plt.grid(True)
    plt.xlabel(symbols[label][var_prev] + r'$_{n-1}$' + ' ' + prev_units)
    plt.ylabel(symbols[label][var_after] + r'$_{n}$' + ' ' + after_units)
    plt.savefig(label + ' ' + var_prev + ' ' + var_after + ' ' + suffix + '.png', dpi=300)
    # plt.show()
    plt.close()


# def previous_division_scatterplots(df, label, var_prev, var_after, suffix=''):
#     prev_array1 = np.concatenate([
#         df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)].sort_values('generation')[var_prev].values[:-1] for dataset in df.dataset.unique() for trap_id in
#         df[df['dataset'] == dataset].trap_ID.unique() for trace
#         in ["A", "B"]
#     ]).flatten()
#     after_array1 = np.concatenate([
#         df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)].sort_values('generation')[var_after].values[1:] for dataset in df.dataset.unique() for trap_id in
#         df[df['dataset'] == dataset].trap_ID.unique() for trace
#         in ["A", "B"]
#     ])
#
#     prev_array = prev_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]
#     after_array = after_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]
#
#     # # Remove outliers
#     # for ar in [prev_array, after_array]:
#     #     before = len(ar)
#     #     ar = ar[np.abs(zscore(ar)) < 3]
#     #     after = len(ar)
#     #     print('percentage that stays after removing all values farther than the fourth standard deviation:', after / before)
#
#     print('done with both')
#
#     pcorr = str(pearsonr(prev_array, after_array)[0])[:5]
#     slope = str(linregress(prev_array, after_array)[0])[:4]
#
#     prev_units = units[var_prev] if label != 'trace_centered' else ''
#     after_units = units[var_after] if label != 'trace_centered' else ''
#
#     sns.set_context('talk')
#
#     fig, ax = plt.subplots(figsize=[3.5, 3.5], tight_layout=True)
#
#     sns.regplot(x=prev_array, y=after_array, data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
#     ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(slope), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
#                 bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
#
#     plt.grid(True)
#     plt.xlabel(symbols[label][var_prev] + r'$_{n-1}$' + ' ' + prev_units)
#     plt.ylabel(symbols[label][var_after] + r'$_{n}$' + ' ' + after_units)
#     plt.savefig(label + ' ' + var_prev + ' ' + var_after + ' ' + suffix + '.png', dpi=300)
#     # plt.show()
#     plt.close()


def main(args):
    def put_all_graphs_into_a_big_grid(df, label, variables=phenotypic_variables, suffix=''):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        unit_symbols = {variable: units[variable] if label != 'trace_centered' else '' for variable in variables}
        
        # Replace the phenotypic variables with their latex counterparts
        df = df[variables].copy().rename(columns=latex_symbols)
        
        # In case we just want to see one scatterplot
        if len(variables) == 2:
    
            sns.set_context('talk')
            
            fig, axes = plt.subplots(figsize=[3.5, 3.5], tight_layout=True)
            
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
            slope, intercept = linregress(relevant[sym2], relevant[sym1])[:2]
            
            axes.scatter(relevant[sym2], relevant[sym1], color='gray', marker='o', alpha=.5)
            axes.plot(np.unique(relevant[sym2]), [intercept + slope * vel for vel in np.unique(relevant[sym2])])
            
            # sns.regplot(x=df[sym2], y=df[sym1], data=df, ax=axes, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
            
            axes.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(str(slope)[:4]), xy=(.5, .79), xycoords=axes.transAxes, fontsize=13, ha='center', va='bottom',
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
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
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
                        
                        ax.scatter(col_array, row_array, color='gray', marker='o', alpha=.5)
                        ax.plot(np.unique(col_array), [intercept + slope * vel for vel in np.unique(col_array)])
                        ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {:.2}$'.format(np.round(slope, 2)), xy=(.5, .72), xycoords=ax.transAxes, fontsize=7, ha='center', va='bottom', color='red',
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
            plt.savefig(args['Figures'] + '/{}/separate/{}/'.format(label, args['data_origin']) + '{}_{}{}.png'.format(variables[0], variables[1], suffix), dpi=300)
        else:
            plt.savefig(args['Figures'] + '/{}/together/'.format(label) + 'regression_plots_{}{}.png'.format(args['data_origin'], suffix), dpi=300)
        # plt.show()
        plt.close()
    
    def heatmap_analogs(df, label, variables=phenotypic_variables, suffix=''):
        
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        
        df = df[variables].copy().rename(columns=latex_symbols)
        
        to_plot = pg.pairwise_corr(df, method='pearson')
        
        # reformat the results into a dataframe we can use for the heatmap
        to_plot_df = pd.DataFrame(columns=to_plot.X.unique(), index=to_plot.Y.unique(), dtype=float)
        repeats = []
        for x in to_plot.X.unique():
            for y in to_plot.Y.unique():
                if x != y and y not in repeats:
                    # Add the results
                    to_plot_df[x].loc[y] = to_plot[(to_plot['X'] == x) & (to_plot['Y'] == y)]['r'].values[0]
                    
            repeats.append(x)
        
        # The mask to show only the lower triangle in the heatmap
        mask = np.ones_like(to_plot_df)
        mask[np.tril_indices_from(mask)] = False
        
        fig, ax = plt.subplots(tight_layout=True, figsize=[7, 5.5])  #  * (2 / 3) on both
        
        sns.heatmap(to_plot_df, annot=True, square=True, vmin=-1, vmax=1, mask=mask, center=0)
        plt.title(args['data_origin'] + ', ' + label)
        # plt.tight_layout()
        plt.savefig(args['Figures'] + '/{}/together/'.format(label) + 'heatmap_{}{}.png'.format(args['data_origin'], suffix), dpi=300)
        # plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scch, label, suffix), dpi=300)
        plt.show()
        plt.close()
    
    # # Create the folder where we will save the figures
    # create_folder('{}/{}'.format(args.figs_location, args.scc))
    # create_folder('{}/{}'.format(args.figs_location, args.scch))
    
    # import the labeled measured bacteria in trace-centered units
    physical_units = pd.read_csv(args['pu'])
    trace_centered = pd.read_csv(args['tc'])
    
    # Calculate this
    time_averages = get_time_averages_df(physical_units, phenotypic_variables)[['lineage_ID', 'max_gen']+phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
    
    for kind, dataframe in zip(['physical_units', 'trace_centered', 'time_averages'], [physical_units, trace_centered, time_averages]):
        print(kind)
        
        # Plot all the variables -- Heatmap
        heatmap_analogs(dataframe, kind, variables=phenotypic_variables, suffix='')

        # Plot all the variables -- Scatter Regression Plot
        put_all_graphs_into_a_big_grid(dataframe, kind, variables=phenotypic_variables, suffix='')

        # Scatter Regression Plots for each variable pair individually
        repeats = []
        for var1 in phenotypic_variables:
            for var2 in phenotypic_variables:
                if var2 in repeats or var2 == var1:
                    continue
                print(var1, var2)
        
                # the separate figures
                put_all_graphs_into_a_big_grid(dataframe, kind, variables=[var1, var2], suffix='')
            repeats.append(var1)
    
    # # Plot all the variables -- Heatmap
    # print('pu')
    # heatmap_analogs(physical_units, 'physical_units', variables=phenotypic_variables, suffix='')
    # print('tc')
    # heatmap_analogs(trace_centered, 'trace_centered', variables=phenotypic_variables, suffix='')
    # print('ta unique')
    # heatmap_analogs(unique_ta, 'time_averages', variables=phenotypic_variables, suffix='')
    #
    # # Plot all the variables -- Scatter Regression Plot
    # print('pu')
    # put_all_graphs_into_a_big_grid(physical_units, 'physical_units', variables=phenotypic_variables[3:], suffix='')
    # print('tc')
    # put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', variables=phenotypic_variables[3:], suffix='')
    # print('unique ta')
    # put_all_graphs_into_a_big_grid(unique_ta, 'time_averages', variables=phenotypic_variables[3:], suffix='')
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
    # parser.add_argument('-kinds_of_correlations', '--kinds_of_correlations', metavar='', nargs="+", help='Calculate pearson and/or variance decomposition correlation?', required=False,
    #                     default=['decomposition', 'pearson'])
    # parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
    #                     default=dict(zip(['main_variables', 'phenotypic_variables'], [main_variables, phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    # Do all the Mother and Sister Machine data
    for data_origin in input_args.dataset_names:
        print(data_origin)
        
        filepath = os.path.dirname(os.path.abspath(__file__))
        
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
            'Figures': filepath+'/Figures',
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
            
        # parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        # parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        # parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
        #                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin)
        # parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
        #                     required=False, default='physical_units.csv')
        # parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
        #                     required=False, default='trace_centered.csv')
        # # parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        # parser.add_argument('-scc', '--scc', metavar='', type=str, help='Where the single cell correlation figures are saved.',
        #                     required=False, default='single_cell_correlations')
        # parser.add_argument('-scch', '--scch', metavar='', type=str, help='Where the single cell correlation heatmap figures are saved.',
        #                     required=False, default='single_cell_correlations_heatmaps')
        # parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
        #                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/'+data_origin)
        #
        # args = parser.parse_args()
        
        main(args)
        
        print('*' * 200)
