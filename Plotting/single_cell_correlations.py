#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, zscore, linregress
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, symbols, units


def previous_division_scatterplots(df, label, var_prev, var_after, suffix=''):
    
    prev_array1 = np.concatenate([
        df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)].sort_values('generation')[var_prev].values[:-1] for dataset in df.dataset.unique() for trap_id in
        df[df['dataset'] == dataset].trap_ID.unique() for trace
        in ["A", "B"]
    ]).flatten()
    after_array1 = np.concatenate([
        df[(df['dataset'] == dataset) & (df['trap_ID'] == trap_id) & (df['trace'] == trace)].sort_values('generation')[var_after].values[1:] for dataset in df.dataset.unique() for trap_id in df[df['dataset'] == dataset].trap_ID.unique() for trace
        in ["A", "B"]
    ])

    prev_array = prev_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]
    after_array = after_array1[(np.abs(zscore(prev_array1)) < 4) & (np.abs(zscore(after_array1)) < 4)]

    # # Remove outliers
    # for ar in [prev_array, after_array]:
    #     before = len(ar)
    #     ar = ar[np.abs(zscore(ar)) < 3]
    #     after = len(ar)
    #     print('percentage that stays after removing all values farther than the fourth standard deviation:', after / before)
    
    print('done with both')
    
    pcorr = str(pearsonr(prev_array, after_array)[0])[:5]
    slope = str(linregress(prev_array, after_array)[0])[:4]

    prev_units = units[var_prev] if label != 'trace_centered' else ''
    after_units = units[var_after] if label != 'trace_centered' else ''
    
    sns.set_context('talk')

    fig, ax = plt.subplots(figsize=[3.5, 3.5], tight_layout=True)

    sns.regplot(x=prev_array, y=after_array, data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
    ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(slope), xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    
    plt.grid(True)
    plt.xlabel(symbols[label][var_prev]+r'$_{n-1}$' + ' ' + prev_units)
    plt.ylabel(symbols[label][var_after] + r'$_{n}$' + ' ' + after_units)
    plt.savefig(label + ' ' + var_prev + ' ' + var_after + ' ' + suffix+'.png', dpi=300)
    # plt.show()
    plt.close()


def main(args):
    def put_all_graphs_into_a_big_grid(df, label, variables=phenotypic_variables, remove_outliers=False, suffix=''):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        unit_symbols = {variable: units[variable] if label != 'trace_centered' else '' for variable in variables}
        
        df = df[variables].copy().rename(columns=latex_symbols)
        if remove_outliers:
            before = len(df)
            df = df[(np.abs(zscore(df)) < 4).all(axis=1)].reset_index(drop=True)
            after = len(df)
            print('percentage that stays after removing all values farther than the fourth standard deviation:', after / before)
        
        if len(variables) == 2:
    
            sns.set_context('talk')
            
            fig, axes = plt.subplots(figsize=[3.5, 3.5], tight_layout=True)

            pcorr = str(pearsonr(df[list(latex_symbols.values())[1]], df[list(latex_symbols.values())[0]])[0])[:4]
            slope = str(linregress(df[list(latex_symbols.values())[1]], df[list(latex_symbols.values())[0]])[0])[:5]
            
            sns.regplot(x=df[list(latex_symbols.values())[1]], y=df[list(latex_symbols.values())[0]], data=df, ax=axes, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
            axes.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(slope), xy=(.5, .72), xycoords=axes.transAxes, fontsize=13, ha='center', va='bottom',
                          color='red',
                          bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
            axes.set_ylabel(list(latex_symbols.values())[0] + ' ' + list(unit_symbols.values())[0])
            axes.set_xlabel(list(latex_symbols.values())[1] + ' ' + list(unit_symbols.values())[1])
            # # plt.tight_layout(pad=.5)
            # plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scc, label, suffix), dpi=300)
            # # plt.show()
            # plt.close()
            
        else:
            
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
                        sns.regplot(x=df[col_var], y=df[row_var], data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
                        ax.annotate(str(pearsonr(df[col_var], df[row_var])[0])[:4], xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
                        # print(df[col_var].quantile([.05, .95]))
                        # exit()
                        # ax.set_xlim(df[col_var].quantile([0, .995]))
                        # ax.set_ylim(df[row_var].quantile([0, .995]))
                        # ax.set_xlim(symbols_bounds[col_var])
                        # ax.set_ylim(symbols_bounds[row_var])
                    
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
        plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scc, label, suffix), dpi=300)
        # plt.show()
        plt.close()
    
    def heatmap_analogs(df, label, variables=phenotypic_variables, suffix=''):
        
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        
        df = df[variables].copy().rename(columns=latex_symbols)
        
        to_plot = df.corr()
        
        mask = np.ones_like(to_plot)
        mask[np.tril_indices_from(mask)] = False
        
        fig, ax = plt.subplots(tight_layout=True, figsize=[7 * (2 / 3), 5.5 * (2 / 3)])
        
        sns.heatmap(to_plot, annot=True, square=True, vmin=-1, vmax=1, mask=mask, center=0)
        # plt.title(label)
        # plt.tight_layout()
        plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scch, label, suffix), dpi=300)
        # plt.show()
        plt.close()
    
    # Create the folder where we will save the figures
    create_folder('{}/{}'.format(args.figs_location, args.scc))
    create_folder('{}/{}'.format(args.figs_location, args.scch))
    
    # import the labeled measured bacteria in trace-centered units
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tc)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    time_averages = pd.read_csv('{}/{}'.format(args.save_folder, args.ta)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    unique_ta = time_averages[['dataset', 'trap_ID', 'trace', 'max_gen', 'fold_growth',
                               'division_ratio', 'added_length', 'generationtime', 'length_birth',
                               'length_final', 'growth_rate']].drop_duplicates().sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])

    previous_division_scatterplots(trace_centered, 'trace_centered', 'division_ratio', 'growth_rate', suffix=' prev and after')
    previous_division_scatterplots(trace_centered, 'trace_centered', 'division_ratio', 'generationtime', suffix=' prev and after')
    previous_division_scatterplots(trace_centered, 'trace_centered', 'division_ratio', 'fold_growth', suffix=' prev and after')
    # exit()
    
    # print()
    
    # heatmap_analogs(trace_centered, 'trace_centered', variables=['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate'], suffix=' main ones')
    # heatmap_analogs(unique_ta, 'unique_ta', variables=['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate'], suffix=' main ones')
    # put_all_graphs_into_a_big_grid(physical_units, 'physical_units', variables=['fold_growth', 'length_birth'], remove_outliers=True, suffix=' phi and x')
    put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', variables=['fold_growth', 'length_birth'], suffix=' phi and x')
    # put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', variables=['fold_growth', 'length_birth'], suffix=' phi and x')
    put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', variables=['generationtime', 'growth_rate'], remove_outliers=True, suffix=' tau and alpha')
    put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', variables=['generationtime', 'growth_rate'], remove_outliers=True, suffix=' tau and alpha')
    # put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', variables=['fold_growth', 'length_birth'], remove_outliers=True, suffix=' phi and ln x')
    
    exit()
    
    print('pu')
    heatmap_analogs(physical_units, 'physical_units', variables=phenotypic_variables)
    print('tc')
    heatmap_analogs(trace_centered, 'trace_centered', variables=phenotypic_variables)
    print('ta')
    heatmap_analogs(time_averages, 'time_averages', variables=phenotypic_variables)
    print('ta unique')
    heatmap_analogs(unique_ta, 'unique_ta', variables=phenotypic_variables)
    
    print('pu')
    put_all_graphs_into_a_big_grid(physical_units, 'physical_units', remove_outliers=True)
    print('tc')
    put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', remove_outliers=True)
    print('ta')
    put_all_graphs_into_a_big_grid(time_averages, 'time_averages', remove_outliers=True)
    print('unique ta')
    put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', remove_outliers=True)

# parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
# parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
# parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units.csv')
# parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages.csv')
# parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered.csv')
# parser.add_argument('-bs', '--bootstraps', metavar='', type=int, help='How many bootstraps per covariance should be done?',
#                     required=False, default=0)
# parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures')
# parser.add_argument('-scc', '--scc', metavar='', type=str, help='Where the single cell correlation figures are saved.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/single_cell_correlations')
# parser.add_argument('-scch', '--scch', metavar='', type=str, help='Where the single cell correlation heatmap figures are saved.',
#                     required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/single_cell_correlations_heatmaps')
#
# args = parser.parse_args()
#
# create_folder(args.save_folder)
# create_folder(args.figs_location)
