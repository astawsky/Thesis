#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, dataset_names, create_folder, phenotypic_variables, cmap, sm_datasets,
    wang_datasets, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets
)
from AnalysisCode.process_rawdata import check_the_division
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

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


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    
""" vd conditioning on trap and lineage """


def vd_with_trap_and_lineage(variables, lin_type, ax):
    # total_df = pd.DataFrame()  # The dataframe with all the experiments in it
    total_df = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/Pooled_SM/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')

    # Check what type of pair lineages we want to look at
    if lin_type == 'NL':
        type_of_lineages = ['NL']
    elif lin_type == 'SL':
        type_of_lineages = ['SL']
    else:
        type_of_lineages = ['SL', 'NL']
        
    total_df = total_df[total_df['dataset'].isin(type_of_lineages)]
    
    print('type of pair lineages we are using:', lin_type)
    output_df = pd.DataFrame()
    
    # What kind of lineage we want to look at, Artificial for noise estimation and Trace for the results
    for kind in ['Trace', 'Artificial']:
        if kind == 'Trace':
            pu = total_df.copy()
            
            # The pooled mean
            pooled_pu_mean = pu[variables].mean()
        else:
            pu = shuffle_info_sm(total_df)  # shuffle the information to get an ergodic set of the same statistics
            
            pooled_pu_mean = pu[variables].mean()  # The pooled mean for the analysis later
        
        # These represent the terms in the equation of the variance decomposition. Here we will store their values inside the sums of the equation of variance for each phenotypic variable seperately.
        delta = pd.DataFrame(columns=variables)
        lin_spec = pd.DataFrame(columns=variables)
        trap = pd.DataFrame(columns=variables)
        
        for trap_id in pu.trap_ID.unique():
            
            trap_ensemble = pu[(pu['trap_ID'] == trap_id)].copy()  # All the bacteria in this experiment
            
            trap_mean = trap_ensemble[variables].mean()  # The average value of all bacteria in this experiment
            
            for trace in ['A', 'B']:
                
                lin = trap_ensemble[(trap_ensemble['trace'] == trace)].copy()  # All the bacteria in this trap
                
                # Add the deviation, for all generations in this lineage, of the mean of this experiment from the mean of the pooled ensemble
                trap = trap.append(lin[variables].count() * ((trap_mean - pooled_pu_mean) ** 2), ignore_index=True)
                
                # Add the deviation, for all generations in this lineage, of the mean of this trap from the mean of this experiment
                lin_spec = lin_spec.append(lin[variables].count() * ((lin[variables].mean() - trap_mean) ** 2), ignore_index=True)
                
                # Add the deviation of all the generation-specific values around the lineage mean
                delta = delta.append(((lin[variables] - lin[variables].mean()) ** 2).sum(), ignore_index=True)
        
        # Calculate the variance for each term with one degree of freedom
        trap_var = trap.sum() / (pu[variables].count() - 1)
        delta_var = delta.sum() / (pu[variables].count() - 1)
        lin_var = lin_spec.sum() / (pu[variables].count() - 1)
        
        # Make sure it is a true decomposition, ie. that the variances of the terms add up to the variance of the pooled ensemble
        assert (np.abs(pu[variables].var() - (trap_var[variables] + delta_var[variables] + lin_var[variables])) < .0000001).all()
        
        # Now save them to the output dataframe
        for variable in variables:
    
            # For every phenotypic variable, add its variance decomposition to the dataframe
            output_df = output_df.append({
                'variable': symbols['physical_units'][variable],  # Latex of the variable (Useful for plotting)
                'Trap+Lin': (trap_var[variable] + lin_var[variable]) / pu[variable].var(),  # Variance explained by the trap and experiment
                'Trap': (trap_var[variable]) / pu[variable].var(),  # Variance explained by the experiment
                'kind': kind
            }, ignore_index=True)

    for spacer in ['', ' ', '  ', '   ', '    ']:
        output_df = output_df.append({
            'variable': spacer,
            'Trap+Lin': 0,
            # This is so we can graph it nicely
            'Trap': 0,
            'kind': 'NA'
        }, ignore_index=True)
    
    real_order = ['', symbols['physical_units']['div_and_fold'], symbols['physical_units']['division_ratio'], symbols['physical_units']['fold_growth'], ' ',
                  symbols['physical_units']['added_length'], symbols['physical_units']['length_birth'], '  ', symbols['physical_units']['generationtime'],
                  symbols['physical_units']['growth_rate'], '   ']
    
    plt.fill_between(output_df.variable.unique(),
                     [output_df[output_df['kind'] == 'Artificial']['Trap+Lin'].mean() for _ in range(len(output_df.variable.unique()))],
                     [0 for _ in range(len(output_df.variable.unique()))], color='lightgrey')
    
    for color, y, label in zip([cmap[0], cmap[1]], ['Trap+Lin', 'Trap'], [r'$\Gamma_{Lin}$', r'$\Gamma_{Env}$']):
        # palette = {"Trace": color, "Artificial": 'red'}
        sns.barplot(x='variable', y=y, data=output_df[output_df['kind'] == 'Trace'], color=color, edgecolor='black', label=label, order=real_order)
    
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, title='')
    
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.ylim([0, .45])


""" vd conditioning on lineages """


def mm_lineage(chosen_datasets, ax):
    total_df = pd.DataFrame()
    for data_origin in chosen_datasets:
        print(data_origin)
        pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(pu)
        pu['experiment'] = data_origin
        
        total_df = total_df.append(pu, ignore_index=True)
    
    output_df = pd.DataFrame(columns=['variable', 'intrinsic', 'environment', 'lineage', 'kind'])
    
    output_df = output_df.append({
        'variable': '',
        # This is so we can graph it nicely
        'Lin': 0,
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
        
        for exp in pu.experiment.unique():
            e_cond = (pu['experiment'] == exp)
            for lin_id in pu[pu['experiment'] == exp].lineage_ID.unique():
                l_cond = (pu['lineage_ID'] == lin_id) & e_cond
                lin = pu[l_cond].copy()
                
                line = line.append(lin[phenotypic_variables].count() * ((lin[phenotypic_variables].mean() - pooled_pu_mean) ** 2), ignore_index=True)
                delta = delta.append(((lin[phenotypic_variables] - lin[phenotypic_variables].mean()) ** 2).sum(), ignore_index=True)
        
        #
        delta_var = delta.sum() / (pu[phenotypic_variables].count() - 1)
        lin_var = line.sum() / (pu[phenotypic_variables].count() - 1)
        
        # Make sure it is a true decomposition
        assert (np.abs(pu[phenotypic_variables].var() - (delta_var[phenotypic_variables] + lin_var[phenotypic_variables])) < .0000001).all()
        
        # Add it to the thing
        for variable in phenotypic_variables:
            output_df = output_df.append({
                'variable': symbols['physical_units'][variable],
                'Lin': (lin_var[variable]) / pu[variable].var(),
                'kind': kind
            }, ignore_index=True)
    
    for kind in ['Trace', 'Artificial']:
        output_df = output_df.append({
            'variable': ' ',
            # This is so we can graph it nicely
            'Lin': 0,
            'kind': kind
        }, ignore_index=True)
        
        output_df = output_df.append({
            'variable': '  ',
            # This is so we can graph it nicely
            'Lin': 0,
            'kind': kind
        }, ignore_index=True)
    
    output_df = output_df.append({
        'variable': '   ',
        # This is so we can graph it nicely
        'Lin': 0,
        'kind': 'Trace'
    }, ignore_index=True)
    
    real_order = ['', symbols['physical_units']['div_and_fold'], symbols['physical_units']['division_ratio'], symbols['physical_units']['fold_growth'], ' ',
                  symbols['physical_units']['added_length'], symbols['physical_units']['length_birth'], '  ', symbols['physical_units']['generationtime'],
                  symbols['physical_units']['growth_rate'], '   ']
    
    conds = (output_df['kind'] == 'Artificial') & (~output_df['variable'].isin(['', ' ', '  ', '   ']))
    
    ax.fill_between(output_df.variable.unique(),
                     [output_df[conds]['Lin'].mean() for _ in range(len(output_df.variable.unique()))],
                     [0 for _ in range(len(output_df.variable.unique()))], color='lightgrey')
    
    # plt.fill_between(real_order,
    #                  [output_df[output_df['kind'] == 'Artificial']['Exp+Lin'].mean() for _ in range(len(real_order))],
    #                  [0 for _ in range(len(real_order))], color='lightgrey')
    
    sns.barplot(x='variable', y='Lin', data=output_df[output_df['kind'] == 'Trace'], color=cmap[0], edgecolor='black', label='Lineage', order=real_order, ax=ax)
    
    handles, labels = ax.get_legend_handles_labels()
    
    # print(handles)
    # print(labels)
    # exit()
    
    ax.legend(handles, [r'$\Gamma_{Trap}$'], title='')
    
    ax.set_xlabel('')
    ax.set_ylabel('Variance Decomposition')
    ax.set_ylim([0, .45])
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


""" vd conditioning on experiments and lineage """


def mm_lineage_experiment(chosen_datasets, ax):
    total_df = pd.DataFrame()
    for data_origin in chosen_datasets:
        print(data_origin)
        pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(pu)
        pu['experiment'] = data_origin
        
        total_df = total_df.append(pu, ignore_index=True)
    
    output_df = pd.DataFrame(columns=['variable', 'intrinsic', 'environment', 'lineage', 'kind'])
    
    output_df = output_df.append({
        'variable': '',
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
    
    for kind in ['Trace', 'Artificial']:
        output_df = output_df.append({
            'variable': ' ',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': kind
        }, ignore_index=True)
        
        output_df = output_df.append({
            'variable': '  ',
            'Exp+Env': 0,
            'Exp+Env+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': kind
        }, ignore_index=True)
    
    output_df = output_df.append({
        'variable': '   ',
        'Exp+Lin': 0,
        # This is so we can graph it nicely
        'Exp': 0,
        'kind': 'Trace'
    }, ignore_index=True)
    
    real_order = ['', symbols['physical_units']['div_and_fold'], symbols['physical_units']['division_ratio'], symbols['physical_units']['fold_growth'], ' ',
                  symbols['physical_units']['added_length'], symbols['physical_units']['length_birth'], '  ', symbols['physical_units']['generationtime'],
                  symbols['physical_units']['growth_rate'], '   ']
    
    conds = (output_df['kind'] == 'Artificial') & (~output_df['variable'].isin(['', ' ', '  ', '   ']))
    
    ax.fill_between(output_df.variable.unique(),
                     [output_df[conds]['Exp+Lin'].mean() for _ in range(len(output_df.variable.unique()))],
                     [0 for _ in range(len(output_df.variable.unique()))], color='lightgrey')
    
    # plt.fill_between(real_order,
    #                  [output_df[output_df['kind'] == 'Artificial']['Exp+Lin'].mean() for _ in range(len(real_order))],
    #                  [0 for _ in range(len(real_order))], color='lightgrey')
    
    for color, y, label in zip([cmap[0], cmap[2]], ['Exp+Lin', 'Exp'], [r'$\Gamma_{Lin}$', r'$\Gamma_{Exp}$']):  # , 'Exp' , 'Experiment'
        # palette = {"Trace": color, "Artificial": 'red'}
        sns.barplot(x='variable', y=y, data=output_df[output_df['kind'] == 'Trace'], color=color, edgecolor='black', label=label, order=real_order, ax=ax)
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='', loc='upper left')
    ax.set_xlabel('')
    ax.set_ylabel('Variance Decomposition')
    ax.set_ylim([0, .45])
    
    
""" vd conditioning on experiment, trap and lineage """


def vd_with_trap_lineage_and_experiments(sm_datasets, variables, lin_type, ax):
    total_df = pd.DataFrame()  # The dataframe with all the experiments in it

    # Check what type of pair lineages we want to look at
    if lin_type == 'NL':
        type_of_lineages = ['NL']
    elif lin_type == 'SL':
        type_of_lineages = ['SL']
    else:
        type_of_lineages = ['SL', 'NL']
    
    # Pool all the experiments and add a category label
    for data_origin in sm_datasets[:-1]:
        print(data_origin)
        
        # pu = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv'
        
        pu = pd.read_csv(os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
        
        pu = pu[pu['dataset'].isin(type_of_lineages)]  # use only the pair lineages we are interested in
        
        if len(total_df) != 0:
            pu['trap_ID'] = pu['trap_ID'] + total_df['trap_ID'].max()
        
        pu['experiment'] = data_origin  # category variable of the experiment
        
        total_df = total_df.append(pu, ignore_index=True)  # Add them to the dataframe with all the experiments
    
    print('type of pair lineages we are using:', lin_type)
    output_df = pd.DataFrame()
    
    # What kind of lineage we want to look at, Artificial for noise estimation and Trace for the results
    for kind in ['Trace', 'Artificial']:
        if kind == 'Trace':
            # pu = total_df[total_df['dataset'].isin(type_of_lineages)].copy()
            pu = total_df.copy()
            
            # The pooled mean
            pooled_pu_mean = pu[variables].mean()
        else:
            pu = shuffle_info_sm(total_df)  # shuffle the information to get an ergodic set of the same statistics
            pu['experiment'] = total_df['experiment'].copy().sort_values().values  # Randomize the experiments as well
            
            pooled_pu_mean = pu[variables].mean()  # The pooled mean for the analysis later
        
        # These represent the terms in the equation of the variance decomposition. Here we will store their values inside the sums of the equation of variance for each phenotypic variable seperately.
        delta = pd.DataFrame(columns=variables)
        diff = pd.DataFrame(columns=variables)
        trap = pd.DataFrame(columns=variables)
        expert = pd.DataFrame(columns=variables)
        
        for exp in pu.experiment.unique():
            
            exp_lins = pu[(pu['experiment'] == exp)].copy()  # All the bacteria in this experiment
            
            exp_mean = exp_lins[variables].mean()  # The average value of all bacteria in this experiment
            
            for trap_id in exp_lins.trap_ID.unique():
                
                trap_lins = exp_lins[(exp_lins['trap_ID'] == trap_id)].copy()  # All the bacteria in this trap
                
                trap_means = trap_lins[variables].mean()  # The average value of all the bacteria in this trap
                
                for trace in ['A', 'B']:
                    
                    lin = trap_lins[(trap_lins['trace'] == trace)].copy()  # All the bacteria in this lineage (channel)
                    
                    # Add the deviation, for all generations in this lineage, of the mean of this experiment from the mean of the pooled ensemble
                    expert = expert.append(lin[variables].count() * ((exp_mean - pooled_pu_mean) ** 2), ignore_index=True)
                    
                    # Add the deviation, for all generations in this lineage, of the mean of this trap from the mean of this experiment
                    trap = trap.append(lin[variables].count() * ((trap_means - exp_mean) ** 2), ignore_index=True)
                    
                    # Add the deviation, for all generations in this lineage, of the mean of this lineage (channel) from the mean of this trap
                    diff = diff.append(lin[variables].count() * ((lin[variables].mean() - trap_means) ** 2), ignore_index=True)
                    
                    # Add the deviation of all the generation-specific values around the lineage mean
                    delta = delta.append(((lin[variables] - lin[variables].mean()) ** 2).sum(), ignore_index=True)
        
        # Calculate the variance for each term with one degree of freedom
        exp_var = expert.sum() / (pu[variables].count() - 1)
        delta_var = delta.sum() / (pu[variables].count() - 1)
        diff_var = diff.sum() / (pu[variables].count() - 1)
        tmean_var = trap.sum() / (pu[variables].count() - 1)
        
        # Make sure it is a true decomposition, ie. that the variances of the terms add up to the variance of the pooled ensemble
        assert (np.abs(pu[variables].var() - (exp_var[variables] + delta_var[variables] + diff_var[variables] + tmean_var[variables])) < .0000001).all()
        
        # Now save them to the output dataframe
        for variable in variables:
    
            # For every phenotypic variable, add its variance decomposition to the dataframe
            output_df = output_df.append({
                'variable': symbols['physical_units'][variable],  # Latex of the variable (Useful for plotting)
                'Exp+Trap': (exp_var[variable] + tmean_var[variable]) / pu[variable].var(),  # Variance explained by the trap and experiment
                'Exp+Trap+Lin': (exp_var[variable] + tmean_var[variable] + diff_var[variable]) / pu[variable].var(),  # Variance explained by the experiment, trap and lineage
                'Exp': (exp_var[variable]) / pu[variable].var(),  # Variance explained by the experiment   if kind == 'Trace' else (exp_var[variable] + tmean_var[variable] + diff_var[variable]) / pu[variable].var()
                'kind': kind
            }, ignore_index=True)

    for spacer in ['', ' ', '  ', '   ', '    ']:
        output_df = output_df.append({
            'variable': spacer,
            'Exp+Trap': 0,
            'Exp+Trap+Lin': 0,
            # This is so we can graph it nicely
            'Exp': 0,
            'kind': 'NA'
        }, ignore_index=True)

    # output_df = output_df.append({
    #     'variable': '',
    #     'Exp+Trap': 0,
    #     'Exp+Trap+Lin': 0,
    #     # This is so we can graph it nicely
    #     'Exp': 0,
    #     'kind': 'NA'
    # }, ignore_index=True)
    #
    # output_df = output_df.append({
    #     'variable': ' ',
    #     'Exp+Trap': 0,
    #     'Exp+Trap+Lin': 0,
    #     # This is so we can graph it nicely
    #     'Exp': 0,
    #     'kind': 'NA'
    # }, ignore_index=True)
    #
    # output_df = output_df.append({
    #     'variable': '  ',
    #     'Exp+Trap': 0,
    #     'Exp+Trap+Lin': 0,
    #     # This is so we can graph it nicely
    #     'Exp': 0,
    #     'kind': 'NA'
    # }, ignore_index=True)
    #
    # output_df = output_df.append({
    #     'variable': '   ',
    #     'Exp+Trap': 0,
    #     'Exp+Trap+Lin': 0,
    #     # This is so we can graph it nicely
    #     'Exp': 0,
    #     'kind': 'NA'
    # }, ignore_index=True)
    #
    # output_df = output_df.append({
    #     'variable': '    ',
    #     'Exp+Trap': 0,
    #     'Exp+Trap+Lin': 0,
    #     # This is so we can graph it nicely
    #     'Exp': 0,
    #     'kind': 'NA'
    # }, ignore_index=True)
    
    real_order = ['', symbols['physical_units']['div_and_fold'], symbols['physical_units']['division_ratio'], symbols['physical_units']['fold_growth'], ' ',
                  symbols['physical_units']['added_length'], symbols['physical_units']['length_birth'], '  ', symbols['physical_units']['generationtime'],
                  symbols['physical_units']['growth_rate'], '   ']
    
    plt.fill_between(output_df.variable.unique(),
                     [output_df[output_df['kind'] == 'Artificial']['Exp+Trap+Lin'].mean() for _ in range(len(output_df.variable.unique()))],
                     [0 for _ in range(len(output_df.variable.unique()))], color='lightgrey')
    
    for color, y, label in zip([cmap[0], cmap[1], cmap[2]], ['Exp+Trap+Lin', 'Exp+Trap', 'Exp'], [r'$\Gamma_{Lin}$', r'$\Gamma_{Trap}$', r'$\Gamma_{Exp}$']):
        # palette = {"Trace": color, "Artificial": 'red'}
        sns.barplot(x='variable', y=y, data=output_df[output_df['kind'] == 'Trace'], color=color, edgecolor='black', label=label, order=real_order)
    
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, title='', loc='upper left')
    
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.ylim([0, .45])
    
    
def cell_cycle_illustration(ax):
    data_origin = 'Pooled_SM'  # We will use cycles from this dataset
    args = {
        'data_origin': data_origin,
        'raw_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/RawData/',
        'processed_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/'
    }
    
    check_the_division(args, lineages=[1], raw_lineages=[], raw_indices=[], pu=[], ax=ax)


scale = 1.5

sns.set_context('paper', font_scale=scale)
sns.set_style("ticks", {'axes.grid': True})

fig, axes = plt.subplots(1, 3, figsize=[6.5 * scale, 3.5 * scale], tight_layout=True)

axes[0].set_title('A', x=-.2, fontsize='xx-large')
axes[1].set_title('B', x=-.3, fontsize='xx-large')
axes[2].set_title('C', x=-.3, fontsize='xx-large')

# axes[0].set_frame_on(False)
# axes[0].set_xticks([])
# axes[0].set_yticks([])
cell_cycle_illustration(axes[0])
axes[0].set_ylabel(r'length $(\mu$m$)$')
axes[0].set_xlabel(r'time $($hr$)$')


mm_lineage(['Lambda_LB', 'Maryam_LongTraces'], ax=axes[1])  # + sm_datasets[:-1], 'MG1655_inLB_LongTraces'
vd_with_trap_and_lineage(phenotypic_variables, lin_type='NL', ax=axes[2])
plt.tight_layout()
plt.show()
plt.close()

exit()

fig, axes = plt.subplots(1, 2, figsize=[6.5 * scale, 3.5 * scale], tight_layout=True)

axes[0].set_title('A', x=-.2, fontsize='xx-large')
axes[1].set_title('B', x=-.05, fontsize='xx-large')

mm_lineage_experiment(['Lambda_LB', 'Maryam_LongTraces'], ax=axes[0])
vd_with_trap_lineage_and_experiments(sm_datasets, phenotypic_variables, lin_type='NL', ax=axes[1])
axes[1].set_yticklabels('')
plt.legend()
plt.show()
plt.close()
