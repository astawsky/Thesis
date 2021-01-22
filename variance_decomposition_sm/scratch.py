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

""" Create all information dataframe where the lineage lengths are kept constant but the cells in the trace itself are randomly sampled from the population without replacement """


def shuffle_info_sm(info):
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


""" decompose wrt traps in a pooled dataset """


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
                pu = shuffle_info_sm(pu)
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
                    'variable': symbols['physical_units'][variable],
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
        plt.savefig(args['save_figs'] + '/' + lin_type, dpi=300)
        # plt.show()
        plt.close()


""" vd conditioning on experiment, trap, and lineage """


def vd_with_trap_and_experiments(sm_datasets, save_fig, variables):
    total_df = pd.DataFrame()
    for data_origin in sm_datasets[:-1]:
        print(data_origin)
        pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/physical_units.csv'
        # pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(pu)
        pu['experiment'] = data_origin
        
        total_df = total_df.append(pu, ignore_index=True)
        
        # print(pu)
        # print(total_df)
        # exit()
    
    # for data_origin in sm_datasets:
    #     pu = os.path.dirname(current_dir) + '/Datasets/' + data_origin + '/ProcessedData/physical_units.csv'
    #     # pu = os.path.dirname(current_dir) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
    #
    
    for type_of_lineages, lin_type in zip([['NL'], ['SL'], ['SL', 'NL']], ['NL', 'SL', 'SL and NL']):
        print('type_of_lineages:', type_of_lineages)
        output_df = pd.DataFrame(columns=['variable', 'intrinsic', 'environment', 'lineage', 'kind'])
        
        output_df = output_df.append({
            'variable': 'before',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': 'Trace'
        }, ignore_index=True)
        
        # Show what datasets it does not have
        if lin_type not in total_df.dataset.unique() and lin_type != 'SL and NL':
            print('does not have {}'.format(lin_type))
            continue
        
        for kind in ['Trace', 'Artificial']:
            if kind == 'Trace':
                # pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
                pu = total_df[total_df['dataset'].isin(type_of_lineages)].copy()
                
                # The pooled mean
                pooled_pu_mean = pu[variables].mean()
            else:
                # pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
                pu = shuffle_info_sm(total_df)
                pu['experiment'] = total_df['experiment'].copy().sort_values().values
                pu = pu[pu['dataset'].isin(type_of_lineages)].copy()
                
                # The pooled mean
                pooled_pu_mean = pu[variables].mean()
            
            #
            delta = pd.DataFrame(columns=variables)
            diff = pd.DataFrame(columns=variables)
            trap = pd.DataFrame(columns=variables)
            expert = pd.DataFrame(columns=variables)
            
            for exp in pu.experiment.unique():
                print('experiment', exp)
                e_cond = (pu['experiment'] == exp)
                exp_lins = pu[e_cond].copy()
                exp_mean = exp_lins[variables].mean()
                for trap_id in pu[pu['experiment'] == exp].trap_ID.unique():
                    t_cond = (pu['trap_ID'] == trap_id) & e_cond
                    trap_lins = pu[t_cond].copy()
                    trap_means = trap_lins[variables].mean()
                    
                    for trace in ['A', 'B']:
                        lin_cond = (trap_lins['trace'] == trace)
                        lin = trap_lins[lin_cond].copy()
                        
                        expert = expert.append(lin[variables].count() * ((exp_mean - pooled_pu_mean) ** 2), ignore_index=True)
                        trap = trap.append(lin[variables].count() * ((trap_means - exp_mean) ** 2), ignore_index=True)
                        diff = diff.append(lin[variables].count() * ((lin[variables].mean() - trap_means) ** 2), ignore_index=True)
                        delta = delta.append(((lin[variables] - lin[variables].mean()) ** 2).sum(), ignore_index=True)
            
            #
            exp_var = expert.sum() / (pu[variables].count() - 1)
            delta_var = delta.sum() / (pu[variables].count() - 1)
            diff_var = diff.sum() / (pu[variables].count() - 1)
            tmean_var = trap.sum() / (pu[variables].count() - 1)
            
            # Make sure it is a true decomposition
            assert (np.abs(pu[variables].var() - (exp_var[variables] + delta_var[variables] + diff_var[variables] + tmean_var[variables])) < .0000001).all()
            
            # Add it to the thing
            for variable in variables:
                
                output_df = output_df.append({
                    'variable': symbols['physical_units'][variable],
                    'Exp+Env': (exp_var[variable] + tmean_var[variable]) / pu[variable].var(),
                    'Exp+Env+Lin': (exp_var[variable] + tmean_var[variable] + diff_var[variable]) / pu[variable].var(),
                    # This is so we can graph it nicely
                    'Exp': (exp_var[variable]) / pu[variable].var() if kind == 'Trace' else (exp_var[variable] + tmean_var[variable] + diff_var[variable]) / pu[variable].var(),
                    'kind': kind
                }, ignore_index=True)
                
                # So we have the spaces in the graph
                if variable == 'fold_growth':
                    output_df = output_df.append({
                        'variable': '',
                        'Exp+Env': 0,
                        'Exp+Env+Lin': 0,
                        # This is so we can graph it nicely
                        'Exp': 0,
                        'kind': kind
                    }, ignore_index=True)
                elif variable == 'length_birth':
                    output_df = output_df.append({
                        'variable': 'between',
                        'Exp+Env': 0,
                        'Exp+Env+Lin': 0,
                        # This is so we can graph it nicely
                        'Exp': 0,
                        'kind': kind
                    }, ignore_index=True)
        
        seaborn_preamble()
        fig, ax = plt.subplots()
        
        output_df = output_df.append({
            'variable': 'after',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': 'Artificial'
        }, ignore_index=True)
        
        plt.fill_between(output_df.variable.unique(),
                         [output_df[output_df['kind'] == 'Artificial']['Exp+Env+Lin'].mean() for _ in range(len(output_df.variable.unique()))],
                         [0 for _ in range(len(output_df.variable.unique()))], color='lightgrey')
        
        for color, y, label in zip([cmap[0], cmap[1], cmap[2]], ['Exp+Env+Lin', 'Exp+Env', 'Exp'], ['Lineage', 'Environment', 'Experiment']):
            # palette = {"Trace": color, "Artificial": 'red'}
            sns.barplot(x='variable', y=y, data=output_df[output_df['kind'] == 'Trace'], color=color, edgecolor='black', label=label)
        
        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles, labels, title='')
        
        plt.title(lin_type)
        plt.xlabel('')
        plt.ylabel('Variance Decomposition')
        plt.ylim([0, .45])
        plt.tight_layout()
        # plt.savefig(save_fig+'/'+lin_type, dpi=300)
        plt.show()
        plt.close()


""" vd conditioning on experiments and lineage """


def vd_with_lineage_and_experiments(chosen_datasets, save_fig, lin_type):
    total_df = pd.DataFrame()
    count = 0
    for data_origin in chosen_datasets:
        print(data_origin)
        pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(pu)
        pu['experiment'] = data_origin
        
        total_df = total_df.append(pu, ignore_index=True)
        
    #     sns.kdeplot(data=pu, x='growth_rate', y='generationtime', label=data_origin, color=cmap[count])
    #     count += 1
    #
    # plt.xlabel(symbols['physical_units']['growth_rate'])
    # plt.ylabel(symbols['physical_units']['generationtime'])
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    # plt.close()

    # exit()
    
    output_df = pd.DataFrame(columns=['variable', 'intrinsic', 'environment', 'lineage', 'kind'])
    
    output_df = output_df.append({
        'variable': 'before',
        'Exp+Env': 0,
        'Exp+Env+Lin': 0,
        # This is so we can graph it nicely
        'Exp': 0,
        'kind': 'Trace'
    }, ignore_index=True)
    
    for kind in ['Trace', 'Artificial']:
        print(kind)
        if kind == 'Trace':
            # pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
            pu = total_df.copy()
            
            # The pooled mean
            pooled_pu_mean = pu[phenotypic_variables].mean()
        else:
            # pu = pd.read_csv(args['pu']).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
            pu = shuffle_info(total_df, mm=True)
            # Shuffle this one manually
            pu['experiment'] = total_df['experiment'].copy().sort_values().values
            pu = pu.copy()
            
            # The pooled mean
            pooled_pu_mean = pu[phenotypic_variables].mean()
        
        #
        delta = pd.DataFrame(columns=phenotypic_variables)
        line = pd.DataFrame(columns=phenotypic_variables)
        expert = pd.DataFrame(columns=phenotypic_variables)
        
        for exp in pu.experiment.unique():
            e_cond = (pu['experiment'] == exp)
            exp_lins = pu[e_cond].copy()
            exp_mean = exp_lins[phenotypic_variables].mean()
            for lin_id in pu[pu['experiment'] == exp].lineage_ID.unique():
                l_cond = (pu['lineage_ID'] == lin_id) & e_cond
                lin = pu[l_cond].copy()
                
                expert = expert.append(lin[phenotypic_variables].count() * ((exp_mean - pooled_pu_mean) ** 2), ignore_index=True)
                line = line.append(lin[phenotypic_variables].count() * ((lin[phenotypic_variables].mean() - exp_mean) ** 2), ignore_index=True)
                delta = delta.append(((lin[phenotypic_variables] - lin[phenotypic_variables].mean()) ** 2).sum(), ignore_index=True)
        
        #
        exp_var = expert.sum() / (pu[phenotypic_variables].count() - 1)
        delta_var = delta.sum() / (pu[phenotypic_variables].count() - 1)
        lin_var = line.sum() / (pu[phenotypic_variables].count() - 1)
        
        # Make sure it is a true decomposition
        assert (np.abs(pu[phenotypic_variables].var() - (exp_var[phenotypic_variables] + delta_var[phenotypic_variables] + lin_var[phenotypic_variables])) < .0000001).all()
        
        # Add it to the thing
        for variable in phenotypic_variables:
            output_df = output_df.append({
                'variable': symbols['physical_units'][variable],
                'Exp+Lin': (exp_var[variable] + lin_var[variable]) / pu[variable].var(),
                # This is so we can graph it nicely
                'Exp': (exp_var[variable]) / pu[variable].var() if kind == 'Trace' else (exp_var[variable] + lin_var[variable]) / pu[variable].var(),
                'kind': kind
            }, ignore_index=True)
    
    seaborn_preamble()
    fig, ax = plt.subplots()
    
    for kind in ['Trace', 'Artificial']:
        output_df = output_df.append({
            'variable': 'between1',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': kind
        }, ignore_index=True)
        
        output_df = output_df.append({
            'variable': 'between2',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': kind
        }, ignore_index=True)
    
    output_df = output_df.append({
        'variable': 'after',
        'Exp+Lin': 0,
        # This is so we can graph it nicely
        'Exp': 0,
        'kind': 'Trace'
    }, ignore_index=True)
    
    real_order = ['before', symbols['physical_units']['div_and_fold'], symbols['physical_units']['division_ratio'], symbols['physical_units']['fold_growth'], 'between1',
                  symbols['physical_units']['added_length'], symbols['physical_units']['length_birth'], 'between2', symbols['physical_units']['generationtime'],
                  symbols['physical_units']['growth_rate'], 'after']
    
    plt.fill_between(real_order,
                     [output_df[output_df['kind'] == 'Artificial']['Exp+Lin'].mean() for _ in range(len(real_order))],
                     [0 for _ in range(len(real_order))], color='lightgrey')
    
    for color, y, label in zip([cmap[0], cmap[1]], ['Exp+Lin', 'Exp'], ['Lineage', 'Experiment']):
        # palette = {"Trace": color, "Artificial": 'red'}
        sns.barplot(x='variable', y=y, data=output_df[output_df['kind'] == 'Trace'], color=color, edgecolor='black', label=label, order=real_order)
    
    handles, labels = ax.get_legend_handles_labels()
    # labels[0] = labels[0] + ': Lin'
    # labels[1] = labels[1] + ': Lin'
    # labels[2] = labels[2] + ': Env'
    # labels[3] = labels[3] + ': Env'
    # labels[4] = labels[4] + ': Exp'
    # labels[5] = labels[5] + ': Exp'
    plt.legend(handles, labels, title='')
    
    plt.title(lin_type)
    plt.xlabel('')
    plt.ylabel('Variance Decomposition')
    # plt.ylim([0, .45])
    plt.tight_layout()
    # plt.savefig(save_fig+'/'+lin_type, dpi=300)
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
    
    kind_of_vd = ['total_length', 'per_gen', 'trap_controlled']
    
    vd_with_lineage_and_experiments(mm_datasets + sm_datasets[:-1], os.path.dirname(os.path.abspath(__file__)), lin_type='mm_datasets')
    exit()
    #
    # vd_with_lineage_and_experiments(cgsc_6300_wang_exps, os.path.dirname(os.path.abspath(__file__)), lin_type='cgsc_6300_wang_exps')
    # exit()
    
    create_folder(os.path.dirname(os.path.abspath(__file__)) + '/Updated')
    
    # vd_with_trap_and_experiments(sm_datasets, os.path.dirname(os.path.abspath(__file__)) + '/Updated',
    #                              variables=['div_and_fold', 'division_ratio', 'fold_growth', 'added_length', 'length_birth', 'generationtime', 'growth_rate'])
    # exit()
    
    # Do all the Mother and Sister Machine data
    for data_origin in sm_datasets[:-1]:  # ['Pooled_SM']:  # wang_datasets:  # input_args.dataset_names:
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
    
    # vd_with_trap(args)
