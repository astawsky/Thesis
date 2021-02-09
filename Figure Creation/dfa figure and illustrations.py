#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os


def plot_dfa_slopes_per_experiment_supplementary():
    se = pd.DataFrame()
    
    for ds in dataset_names:
        if ds == 'Pooled_SM':
            continue
        print(ds)
        scaling_exponents = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/DFA/Scaling Exponents/{ds}/scaling_exponents.csv')
        scaling_exponents['experiment'] = ds
        
        se = se.append(scaling_exponents, ignore_index=True)
    
    se = se.replace({'$f_n e^{\phi_{n}}$': r'$r$'})
    
    condition = (se['dataset'] == 'Trace') & (se['kind'] == 'dfa (short)')
    # condition2 = (se['dataset'] == 'Shuffled') & (se['kind'] == 'dfa (short)')
    
    sns.set_context('paper', font_scale=1.5)
    sns.set_style("ticks", {'axes.grid': True})
    fig, ax = plt.subplots(figsize=[7, 7], tight_layout=True)
    plt.axhline(0.5, ls='-', c='k')
    sns.pointplot(data=se[condition], x='variable', y='slope', hue='experiment', join=False, dodge=True, palette=cmap, ci="sd")
    # sns.boxplot(data=se, x='variable', y='slope', hue='dataset', showfliers=False, palette=cmap)
    # sns.pointplot(data=se, x='variable', y='slope', hue='dataset', join=False, dodge=True, palette=cmap, ci="sd", capsize=.1)
    # sns.violinplot(data=se, x='variable', y='slope', hue='dataset', join=False, dodge=True, color=cmap[0])
    # plt.legend(title='')
    plt.ylabel(r'$\gamma$')
    plt.xlabel('')
    ax.get_legend().remove()
    plt.ylim([0, 1.2])
    plt.savefig('Figures/supplemental, dfa slopes per experiment.png', dpi=300)
    # plt.show()
    plt.close()
    
    
def plot_dfa_cumsum_illustration(ax):
    min_ws = 5
    max_ws = None
    window_size_steps = 2
    steps_between_windows = 3
    c = cmap[3]
    
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    
    lin_id = 15
    rel = pu[pu['lineage_ID'] == lin_id].copy().sort_values('generation')
    
    r = rel['length_birth']
    r = (r - r.mean()).cumsum().values
    
    ws = 5
    
    starting = np.arange(0, len(r) - ws, 3)
    ending = starting + ws - 1
    
    ax.axhline(0, c='k')
    
    for start, end in zip(starting, ending):
        ax.axvline(start + 1, color=cmap[5])
        ax.axvline(end + 1, color=cmap[6])

    ax.plot(np.arange(1, len(r) + 1), r, label=r'$Z$')

    for start, end in zip(starting, ending):
        slope, intercept = linregress(np.arange(start + 1, end + 2), r[start:end + 1])[:2]
        ax.plot(np.arange(start + 1, end + 2), intercept + slope * np.arange(start + 1, end + 2), color=c, ls='--', label='' if start != starting[-1] else r'$\hat{Z}$')
    
    ax.set_xlim([1, len(r)])
    ax.set_xlabel('Generation')
    ax.legend(title='')
    # plt.savefig('Figures/', dpi=300)
    # plt.show()
    # plt.close()


def recreate_loglog_plots(loglog_plot_recreation, variable, ax):
    print(variable)
    
    # For every kind of loglog plot
    relevant = loglog_plot_recreation[(loglog_plot_recreation['variable'] == variable) & (loglog_plot_recreation['dataset'].isin(['Trace', 'Shuffled']))].copy()
    
    set_trace_legend, set_shuffled_legend = False, False
    
    wss_array = {'Trace': np.array([]), 'Shuffled': np.array([])}
    ytr_array = wss_array.copy()
    
    # For each lineage separately
    for ind in relevant.index:
        # Convert the window sizes and mean squared error per windows size from the dataframe from a string to an array of integers or floats
        wss = relevant['window_sizes'].loc[ind].strip('][').split(' ')
        wss = [int(r.split('\n')[0]) for r in wss if r != '']
        ytr = relevant['y_to_regress'].loc[ind].strip('][').split(', ')
        ytr = [float(r) for r in ytr if r != '']
        
        # If it is a trace lineage plot it and include it in the legend
        if relevant.loc[ind]['dataset'] == 'Trace':
            # The constant shuffled color
            color = cmap[0]
            
            # Add the analysis curve to the array of its dataset for the population regression
            wss_array['Trace'] = np.append(wss_array['Trace'], wss)
            ytr_array['Trace'] = np.append(ytr_array['Trace'], ytr)
            
            # Include it in the legend but not more than once
            if set_trace_legend:
                ax.plot(wss, ytr, color=color, marker='')
            else:
                ax.plot(wss, ytr, color=color, marker='')  # , label='Trace'
                set_trace_legend = True
        # If it is a shuffled lineage plot it and include it in the legend
        elif relevant.loc[ind]['dataset'] == 'Shuffled':
            # The constant shuffled color
            color = cmap[1]
            
            # Add the analysis curve to the array of its dataset for the population regression
            wss_array['Shuffled'] = np.append(wss_array['Shuffled'], wss)
            ytr_array['Shuffled'] = np.append(ytr_array['Shuffled'], ytr)
            
            # Include it in the legend but not more than once
            if set_shuffled_legend:
                ax.plot(wss, ytr, color=color, marker='')
            else:
                ax.plot(wss, ytr, color=color, marker='')  # , label='Shuffled'
                set_shuffled_legend = True
        # We do not want to see the white noise
        else:
            continue
    
    # Get the linear regression of all the trajectories pooled for each dataset
    slope_trace, intercept_trace, _, _, std_err_trace = linregress(np.log(np.array(wss_array['Trace']).flatten()), np.log(np.array(ytr_array['Trace']).flatten()))
    slope_art, intercept_art, _, _, std_err_art = linregress(np.log(np.array(wss_array['Shuffled']).flatten()), np.log(np.array(ytr_array['Shuffled']).flatten()))
    
    # Plot the best fit line and its parameters
    ax.plot(np.unique(wss_array['Trace']), [np.exp(intercept_trace) * (l ** slope_trace) for l in np.unique(wss_array['Trace'])], ls='--', color='blue', linewidth=3,
             label='Trace', zorder=100)  # label=r'$' + str(np.round(intercept_trace, 2)) + r'\, k^{' + str(np.round(slope_trace, 2)) + r'\pm' + str(np.round(std_err_trace, 3)) + r'}$'
    ax.plot(np.unique(wss_array['Shuffled']), [np.exp(intercept_art) * (l ** slope_art) for l in np.unique(wss_array['Shuffled'])], ls='--', color='red', linewidth=3,
             label='Shuffled', zorder=100)  # label=r'$' + str(np.round(intercept_art, 2)) + r'\, k^{' + str(np.round(slope_art, 2)) + r'\pm' + str(np.round(std_err_art, 3)) + r'}$'
    
    # ax.title(variable)
    ax.set_ylabel(r'$F(k)$')
    ax.set_xlabel('k (window size)')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend(title='', loc='upper left')


def plot_dfa_slopes(ax):
    se = pd.DataFrame()
    
    for ds in dataset_names[:-1]:
        if ds == 'Pooled_SM':
            continue
        print(ds)
        scaling_exponents = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/DFA/Scaling Exponents/{ds}/scaling_exponents.csv')
        scaling_exponents['experiment'] = ds
        
        se = se.append(scaling_exponents, ignore_index=True)
    
    se = se.replace({'$f_n e^{\phi_{n}}$': r'$r$'})
    
    # condition = (se['dataset'] == 'Trace') & (se['kind'] == 'dfa (short)')
    # condition2 = (se['dataset'] == 'Shuffled') & (se['kind'] == 'dfa (short)')
    
    plt.axhline(0.5, ls='-', c='k')
    # sns.pointplot(data=se[condition], x='variable', y='slope', hue='experiment', join=False, dodge=True, palette=cmap, ci="sd")
    # sns.boxplot(data=se, x='variable', y='slope', hue='dataset', showfliers=False, palette=cmap)
    
    sns.pointplot(data=se, x='variable', y='slope', hue='dataset', join=False, dodge=True, palette=cmap, capsize=.1, ax=ax, ci="sd", zorder=100)
    
    # for count, datas in enumerate(se.dataset.unique()):
    #     c = cmap[count]
    #     y = [se[(se['dataset'] == datas) & (se['variable'] == v)]['slope'].mean() for v in se.variable.unique()]
    #     yerr = [se[(se['dataset'] == datas) & (se['variable'] == v)]['slope'].std() for v in se.variable.unique()]
    #     ax.errorbar(x=se.variable.unique(), y=y, yerr=yerr, capsize=.1, c=c, label=datas)
    
    # sns.violinplot(data=se, x='variable', y='slope', hue='dataset', join=False, dodge=True, color=cmap[0])
    # plt.legend(title='')
    ax.set_ylabel(r'$\gamma$')
    ax.set_xlabel('')
    ax.legend(title='')
    # ax.get_legend().remove()
    plt.ylim([0, 1.2])
    

# plot_dfa_slopes_per_experiment_supplementary()

recreation = pd.read_csv(r'/Users/alestawsky/PycharmProjects/Thesis/DFA/LogLog Recreation/MG1655_inLB_LongTraces/loglog_scaling_recreation.csv')

sns.set_context('paper', font_scale=1)
sns.set_style("ticks", {'axes.grid': True})
fig, axes = plt.subplots(1, 3, figsize=[6.5, 3.5], tight_layout=True)

axes[0].set_title('A', x=-.1, fontsize='xx-large')
axes[1].set_title('B', x=-.12, fontsize='xx-large')
axes[2].set_title('C', x=-.18, fontsize='xx-large')

plot_dfa_cumsum_illustration(axes[0])
recreate_loglog_plots(recreation[recreation['kind'] == 'dfa (short)'], '$x_0$', axes[1])
plot_dfa_slopes(axes[2])
plt.savefig('Figures/dfa figure.png', dpi=500)
# plt.show()
plt.close()
