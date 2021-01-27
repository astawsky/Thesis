#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import linregress, pearsonr, spearmanr
import pingouin as pg
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from AnalysisCode.global_variables import phenotypic_variables, create_folder, symbols, units, dataset_names, get_time_averages_df, sm_datasets, wang_datasets, tanouchi_datasets, seaborn_preamble


def previous_division_scatterplots(df, label, var_prev, var_after, suffix=''):
    prev_array = np.concatenate([
        df[(df['lineage_ID'] == lineage_id)].sort_values('generation')[var_prev].values[:-1] for lineage_id in df.lineage_ID.unique()]).flatten()
    after_array = np.concatenate([
        df[(df['lineage_ID'] == lineage_id)].sort_values('generation')[var_after].values[1:] for lineage_id in df.lineage_ID.unique()]).flatten()
    
    pcorr = str(pg.corr(prev_array, after_array, method='spearman')['r'].loc['spearman'])[:4]
    slope = str(linregress(prev_array, after_array)[0])[:4]
    
    prev_units = units[var_prev] if label != 'trace_centered' else ''
    after_units = units[var_after] if label != 'trace_centered' else ''
    
    sns.set_context('talk')
    
    fig, ax = plt.subplots(figsize=[5, 5], tight_layout=True)
    
    sns.regplot(x=prev_array, y=after_array, data=df, ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
    ax.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(slope), xy=(.5, .92), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    
    plt.grid(True)
    plt.xlabel(symbols[label][var_prev] + r'$_{n-1}$' + ' ' + prev_units)
    plt.ylabel(symbols[label][var_after] + r'$_{n}$' + ' ' + after_units)
    plt.savefig(label + ' ' + var_prev + ' ' + var_after + ' ' + suffix + '.png', dpi=300)
    plt.show()
    plt.close()


def main(df, label, var1, var2):
    
    variables = [var1, var2]
    
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
    df = df[variables + ['experiment']].copy().rename(columns=latex_symbols)
    
    # stylistic reasons
    sns.set_context('paper', font_scale=1.2)
    sns.set_style("ticks", {'axes.grid': True})
    
    fig, axes = plt.subplots(figsize=[5, 5], tight_layout=True)
    
    sym1 = list(latex_symbols.values())[0]
    sym2 = list(latex_symbols.values())[1]
    
    relevant = df[[sym2, sym1, 'experiment']].dropna()
    
    pcorr = str(pearsonr(relevant[sym2].values, relevant[sym1].values)[0])[:4]
    slope, intercept, r_value, _, std_err = linregress(relevant[sym2], relevant[sym1])
    
    # sns.regplot(data=relevant, x=sym2, y=sym1, line_kws={'color': 'black', 'ls': '--', 'lw': 2})
    sns.scatterplot(data=relevant, x=sym2, y=sym1, hue='experiment')
    plt.legend(title='')
    
    # sns.regplot(x=df[sym2], y=df[sym1], data=df, ax=axes, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
    
    axes.annotate(r'$\rho = $' + pcorr + '\n' + r'$\beta = {}$'.format(str(slope)[:4]), xy=(.5, .92), xycoords=axes.transAxes, fontsize=13, ha='center', va='bottom',
                  color='red',
                  bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
    axes.set_ylabel(sym1 + ' ' + list(unit_symbols.values())[0])
    axes.set_xlabel(sym2 + ' ' + list(unit_symbols.values())[1])
    
    plt.savefig(r'/Users/alestawsky/PycharmProjects/Thesis/sm_ta_correlations_figure/experiment_centered_gr.png', dpi=300)
    # plt.show()
    plt.close()


if __name__ == '__main__':
    import argparse
    import os
    
    # The variables we want to plot
    main_variables = ['fold_growth', 'division_ratio', 'generationtime', 'length_birth', 'growth_rate']

    total_pu = pd.DataFrame()
    total_tc = pd.DataFrame()
    total_ta = pd.DataFrame()
    for data_origin in sm_datasets[:-1]:
        print(data_origin)
        pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/physical_units.csv')
        tc = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/trace_centered.csv')
        
        ta = get_time_averages_df(pu, phenotypic_variables)[['lineage_ID', 'max_gen'] + phenotypic_variables].drop_duplicates().sort_values(['lineage_ID'])
        
        pool = pu[phenotypic_variables].mean()
        
        pu['experiment'] = data_origin
        pu[phenotypic_variables] = pu[phenotypic_variables] - pool
        total_pu = total_pu.append(pu, ignore_index=True)
        tc['experiment'] = data_origin
        tc[phenotypic_variables] = tc[phenotypic_variables]
        total_tc = total_tc.append(tc, ignore_index=True)
        ta['experiment'] = data_origin
        ta[phenotypic_variables] = ta[phenotypic_variables] - pool
        total_ta = total_ta.append(ta, ignore_index=True)
    
    main(total_ta, 'time_averages', 'growth_rate', 'growth_rate')
