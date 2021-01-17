#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress
from itertools import combinations

""" Create all information dataframe where the lineage lengths are kept constant but the cells in the trace itself are randomly sampled from the population without replacement """


def shuffle_info(info):
    # Give it a name, contains S, NL
    new_info = pd.DataFrame(columns=info.columns)
    
    # what is the trace length of each trace? This is the only thing that stays the same
    sizes = {'{} {} {}'.format(dataset, trap_ID, trace): len(info[(info['trap_ID'] == trap_ID) & (info['trace'] == trace) & (info['dataset'] == dataset)]) for dataset in np.unique(info['dataset'])
             for trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']) for trace in ['A', 'B']}
    
    lineage_id = 0
    for dataset in np.unique(info['dataset']):
        for trap_ID in np.unique(info[info['dataset'] == dataset]['trap_ID']):
            for trace in ['A', 'B']:
                # trace length
                size = sizes['{} {} {}'.format(dataset, trap_ID, trace)]
                
                # sample from the old dataframe
                samples = info[info['dataset'] == dataset].sample(replace=False, n=size)
                
                # drop what we sampled
                info = info.drop(index=samples.index)
                
                # add some correct labels even though they don't matter so that the dataframe structure is still intact
                samples['dataset'] = dataset
                samples['trap_ID'] = trap_ID
                samples['trace'] = trace
                samples['lineage_ID'] = lineage_id
                samples['generation'] = np.arange(size)
                
                # add them to the new, shuffled dataframe
                new_info = new_info.append(samples, ignore_index=True)
                
                # Go on to the next ID
                lineage_id += 1

    return new_info


def change_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


""" vd conditioning on trap as well as lineage """


def vd_with_trap(args):
    
    for type_of_lineages, lin_type in zip([['SL'], ['NL'], ['SL', 'NL']], ['SL', 'NL', 'SL and NL']):
        output_df = pd.DataFrame(columns=['variable', 'intrinsic', 'environment', 'lineage', 'kind'])
        
        if lin_type not in pd.read_csv(args['pu']).dataset.unique() and lin_type != 'SL and NL':
            print('does not have {}'.format(lin_type))
            continue
        
        for kind in ['Trace', 'Artificial']:
            if kind == 'Trace':
                pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
                pu = pu[pu['dataset'].isin(type_of_lineages)].copy()
        
                # The pooled mean
                pooled_pu_mean = pu[phenotypic_variables].mean()
            else:
                pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
                pu = shuffle_info(pu)
                pu = pu[pu['dataset'].isin(type_of_lineages)].copy()
        
                # The pooled mean
                pooled_pu_mean = pu[phenotypic_variables].mean()
    
            #
            delta = pd.DataFrame(columns=phenotypic_variables)
            diff = pd.DataFrame(columns=phenotypic_variables)
            trap = pd.DataFrame(columns=phenotypic_variables)
            
            for trap_id in pu.trap_ID.unique():
                t_cond = (pu['trap_ID'] == trap_id)
                trap_lins = pu[t_cond].copy()
                trap_means = trap_lins[phenotypic_variables].mean()
                
                for trace in ['A', 'B']:
                    lin_cond = (trap_lins['trace'] == trace)
                    lin = trap_lins[lin_cond].copy()
                    
                    trap = trap.append(lin[phenotypic_variables].count() * ((trap_means - pooled_pu_mean) ** 2), ignore_index=True)
                    diff = diff.append(lin[phenotypic_variables].count() * ((trap_means - lin[phenotypic_variables].mean()) ** 2), ignore_index=True)
                    delta = delta.append(((lin[phenotypic_variables] - lin[phenotypic_variables].mean()) ** 2).sum(), ignore_index=True)
                    
            #
            delta_var = delta.sum() / (pu[phenotypic_variables].count() - 1)
            diff_var = diff.sum() / (pu[phenotypic_variables].count() - 1)
            tmean_var = trap.sum() / (pu[phenotypic_variables].count() - 1)
            
            # Add it to the thing
            for variable in phenotypic_variables:
                output_df = output_df.append({
                    'variable': symbols['time_averages'][variable],
                    'environment': (tmean_var[variable]) / pu[variable].var(),
                    'environment+lineage': (tmean_var[variable] + diff_var[variable]) / pu[variable].var(),
                    'kind': kind
                }, ignore_index=True)
        
        seaborn_preamble()
        fig, ax = plt.subplots()
        for color, y in zip([cmap[0], cmap[1]], ['environment+lineage', 'environment']):
            palette = {"Trace": color, "Artificial": change_color(color)}
            sns.barplot(x='variable', y=y, data=output_df, hue='kind', palette=palette, edgecolor='black')

        handles, labels = ax.get_legend_handles_labels()
        labels[0] = labels[0] + ': lineage'
        labels[1] = labels[1] + ': lineage'
        labels[2] = labels[2] + ': environment'
        labels[3] = labels[3] + ': environment'
        plt.legend(handles, labels, title='')
        
        plt.title(lin_type)
        plt.xlabel('')
        plt.ylabel('Variance Decomposition')
        plt.ylim([0, .45])
        plt.tight_layout()
        plt.savefig(args['save_figs']+'/'+lin_type, dpi=300)
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
    parser.add_argument('-kinds_of_correlations', '--kinds_of_correlations', metavar='', nargs="+", help='Calculate pearson and/or variance decomposition correlation?', required=False,
                        default=['decomposition', 'pearson'])
    parser.add_argument('-variable_mapping', '--variable_mapping', metavar='', nargs="+", help='Calculate for what variables in the figure?', required=False,
                        default=dict(zip(['phenotypic_variables'], [phenotypic_variables])))
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    kind_of_vd = ['total_length', 'per_gen', 'trap_controlled']
    
    # Do all the Mother and Sister Machine data
    for data_origin in sm_datasets:#['Pooled_SM']:  # wang_datasets:  # input_args.dataset_names:
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
            'tc': processed_data + 'z_score_under_3/trace_centered_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'trace_centered.csv'
        }
        
        vd_with_trap(args)
