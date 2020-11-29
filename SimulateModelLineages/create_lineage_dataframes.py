#!/usr/bin/env bash

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder
import argparse
import os
import time
from itertools import combinations
pd.options.mode.chained_assignment = None  # default='warn', FALSE POSITIVE WARNING

first_time = time.time()
start_time = first_time

parser = argparse.ArgumentParser(description='Process Lineage Data.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                    required=False, default='physical_units.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                    required=False, default='trace_centered.csv')


def conditional_normal_parameters(a, mean_1, mean_2, sigma_11, sigma_22, sigma_12):
    sigma_21 = sigma_12.T
    
    if len(sigma_22) == 1:
        inverse_matrix = 1 / sigma_22.values[0][0]
    else:
        inverse_matrix = pd.DataFrame(np.linalg.pinv(sigma_22.values), sigma_22.columns, sigma_22.index)
    
    if (len(a) == 1) and (len(sigma_22) == 1):
        
        new_mean = mean_1.values[0] + sigma_12.values[0][0] * inverse_matrix * (a.values[0] - mean_2.values[0])
        
        new_cov = sigma_11.values[0][0] - sigma_12.values[0][0] * inverse_matrix * sigma_21.values[0][0]
        
        freshly_sampled = np.random.normal(new_mean, np.sqrt(new_cov), 1)
    else:
        
        new_mean = mean_1 + sigma_12 @ inverse_matrix @ (a - mean_2)  # np.linalg.inv(pd.DataFrame(np.linalg.pinv(sigma_22.values), sigma_22.columns, sigma_22.index))
        
        # print(new_mean)
        
        new_cov = sigma_11 - sigma_12 @ inverse_matrix @ sigma_21
        # print(new_cov)
        
        freshly_sampled = np.random.multivariate_normal(new_mean, new_cov, 1)[0]
    
    return freshly_sampled


def get_inter_covariance(trace_centered):
    # Where we will put the rows to calculate the covariance matrix
    old = []
    young = []
    
    for dataset in np.unique(trace_centered['dataset']):
        for trap_id in np.unique(trace_centered[trace_centered['dataset'] == dataset]['trap_ID']):
            for trace in ['A', 'B']:
                tc_lineage = trace_centered[(trace_centered['dataset'] == dataset) & (trace_centered['trap_ID'] == trap_id) & (trace_centered['trace'] == trace)].copy()
                for gen in np.unique(tc_lineage['generation'])[:-1]:
                    young.append(tc_lineage[tc_lineage['generation'] == gen])
                    old.append(tc_lineage[tc_lineage['generation'] == gen + 1])
    
    previous_variables = ['previous, ' + variable for variable in phenotypic_variables]
    
    inter_covariance = pd.DataFrame(index=phenotypic_variables, columns=previous_variables)
    
    for var1, v1 in zip(phenotypic_variables, previous_variables):
        print(var1)
        for var2 in phenotypic_variables:
            inter_covariance[v1].loc[var2] = float(np.cov([float(oldie[var1]) for oldie in old], [float(youngie[var2]) for youngie in young])[0, 1])
    
    return inter_covariance


def based_on_target(targets, dep, pu_lineage, sg_cov, inter_covariance):
    sampled_lineage = pu_lineage[pu_lineage['generation'] == 0].copy()
    
    cov = pd.DataFrame(columns=dep, index=targets)
    mean2 = pd.Series(index=dep)
    a = pd.Series(index=dep)
    
    sigma_22 = pd.DataFrame(index=dep, columns=dep, dtype=float)
    
    for variable in dep:
        # A same gen variable
        if len(variable.split(', ')) == 1:
            cov[variable].loc[targets] = sg_cov[variable].loc[targets]
            mean2[variable] = pu_lineage[variable].mean()
            a[variable] = pu_lineage[pu_lineage['generation'] == 1][variable].values
            
            for variable_col in dep:
                if len(variable_col.split(', ')) == 1:
                    sigma_22[variable_col].loc[variable] = sg_cov[variable_col].loc[variable]
                else:
                    sigma_22[variable_col].loc[variable] = sg_cov[variable_col.split(', ')[1]].loc[variable]
        else:
            cov[variable].loc[targets] = inter_covariance[variable].loc[targets]
            mean2[variable] = pu_lineage[variable.split(', ')[1]].mean()
            a[variable] = pu_lineage[pu_lineage['generation'] == 0][variable.split(', ')[1]].values
            
            for variable_col in dep:
                if len(variable_col.split(', ')) == 1:
                    sigma_22[variable_col].loc[variable] = sg_cov[variable_col].loc[variable.split(', ')[1]]
                else:
                    sigma_22[variable_col].loc[variable] = sg_cov[variable_col.split(', ')[1]].loc[variable.split(', ')[1]]
    
    for gen in np.sort(pu_lineage['generation'])[1:]:  # From the previous generation
        
        to_add = {col: [] for col in pu_lineage.columns}
        
        gen = int(gen)
        
        freshly_sampled = conditional_normal_parameters(a, pu_lineage[targets].mean(), mean2, sg_cov[targets].loc[targets], sigma_22, cov)
        
        for target, sample in zip(targets, freshly_sampled):
            to_add[target].append(sample)
        
        to_add['generation'].append(gen)
        
        if targets == ['generationtime']:
            # We sample this randomly, which is almost as if it's constant
            gr = np.random.normal(pu_lineage['growth_rate'].mean(), pu_lineage['growth_rate'].std(), 1)[0]
            dr = np.random.normal(pu_lineage['division_ratio'].mean(), pu_lineage['division_ratio'].std(), 1)[0]
            lb = float(sampled_lineage[sampled_lineage['generation'] == (gen - 1)]['division_ratio'] * sampled_lineage[sampled_lineage['generation'] == (gen - 1)]['length_final'])
            to_add['growth_rate'].append(gr)
            to_add['division_ratio'].append(dr)
            # We get this from the mapping
            to_add['length_birth'].append(lb)
            to_add['length_final'].append(float(lb * np.exp(gr * to_add['generationtime'][-1])))
            to_add['added_length'].append(float(to_add['length_final'][-1] - to_add['length_birth'][-1]))
            to_add['fold_growth'].append(float(gr * to_add['generationtime'][-1]))
        elif targets == ['growth_rate']:
            # We sample this randomly, which is almost as if it's constant
            gt = np.random.normal(pu_lineage['generationtime'].mean(), pu_lineage['generationtime'].std(), 1)[0]
            dr = np.random.normal(pu_lineage['division_ratio'].mean(), pu_lineage['division_ratio'].std(), 1)[0]
            lb = float(sampled_lineage[sampled_lineage['generation'] == (gen - 1)]['division_ratio'] * sampled_lineage[sampled_lineage['generation'] == (gen - 1)]['length_final'])
            to_add['generationtime'].append(gt)
            to_add['division_ratio'].append(dr)
            # We get this from the mapping
            to_add['length_birth'].append(lb)
            to_add['length_final'].append(float(lb * np.exp(gt * to_add['growth_rate'][-1])))
            to_add['added_length'].append(float(to_add['length_final'][-1] - to_add['length_birth'][-1]))
            to_add['fold_growth'].append(float(gt * to_add['growth_rate'][-1]))
        else:  # meaning we have both generationtime and growth rate
            dr = np.random.normal(pu_lineage['division_ratio'].mean(), pu_lineage['division_ratio'].std(), 1)[0]
            lb = sampled_lineage[sampled_lineage['generation'] == gen - 1]['division_ratio'] * sampled_lineage[sampled_lineage['generation'] == gen - 1]['length_final']
            to_add['growth_rate'].append(to_add['generationtime'][-1])
            to_add['division_ratio'].append(dr)
            # We get this from the mapping
            to_add['length_birth'].append(lb)
            to_add['length_final'].append(lb * np.exp(to_add['generationtime'][-1] * to_add['growth_rate'][-1]))
            to_add['added_length'].append(to_add['length_final'][-1] - to_add['length_birth'][-1])
            to_add['fold_growth'].append(to_add['generationtime'][-1] * to_add['growth_rate'][-1])
        
        # The categorical variables
        to_add['dataset'].append(np.unique(pu_lineage['dataset'])[0])
        to_add['trap_ID'].append(np.unique(pu_lineage['trap_ID'])[0])
        to_add['trace'].append(np.unique(pu_lineage['trace'])[0])
        
        # update the conditions
        for variable in dep:
            if len(variable.split(', ')) == 1:
                a[variable] = to_add[variable][-1]
            else:
                a[variable] = to_add[variable.split(', ')[1]][-1]
        
        # add them to the lineage we want to create
        sampled_lineage = sampled_lineage.append(pd.DataFrame(to_add), ignore_index=True)
    
    # print('\n' * 5)
    # print(sampled_lineage['length_birth'])
    # print(sampled_lineage['generationtime'])
    # print(sampled_lineage['growth_rate'])
    # print(sampled_lineage['length_final'])
    # print(sampled_lineage['division_ratio'])
    # print(sampled_lineage['fold_growth'])
    # exit()
    
    return sampled_lineage


def main(args):
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tc))
    
    # In case we have not calculated it: young rows, old columns
    try:
        inter_covariance = pd.read_csv('{}/{}'.format(args.save_folder, args.inter_cov), index_col=0)
    except:
        inter_covariance = get_inter_covariance(trace_centered)
        inter_covariance.to_csv('{}/{}'.format(args.save_folder, args.inter_cov), index=True)
    
    possible_dependencies = ['previous, ' + variable for variable in ['generationtime', 'growth_rate', 'division_ratio']] + ['length_birth']
    
    actual_dependencies = []
    for n in np.arange(1, len(possible_dependencies)):
        actual_dependencies += [list(val) for val in list(combinations(possible_dependencies, n))]
    
    sg_cov = trace_centered[['generationtime', 'growth_rate', 'length_birth', 'division_ratio']].cov()
    
    # for targets in [['generationtime'], ['growth_rate'], ['generationtime', 'growth_rate']]:
    #     for dep in actual_dependencies:
    #         # print('targets:', targets, 'dependencies:', dep)
    #         filename = 'target ' + ', '.join(targets) + ' and dependent on ' + ', '.join(
    #             ['same_' + variable.split(', ')[0] if len(variable.split(', ')) == 1 else variable.split(', ')[1] for variable in dep])
    #         print(filename)
    #
    # exit()
    
    for targets in [['generationtime', 'growth_rate']]: # ['generationtime'], ['growth_rate'],
        for dep in actual_dependencies:
            # targets = ['generationtime']
            # dep = ['previous, generationtime', 'previous, growth_rate']

            filename = 'target ' + ', '.join(targets) + ' and dependent on ' + ', '.join(
                ['same_' + variable.split(', ')[0] if len(variable.split(', ')) == 1 else variable.split(', ')[1] for variable in dep])
            
            print('targets:', targets, 'dependencies:', dep)
            
            # These are the dataframes that we will output
            to_output_physical_units = pd.DataFrame(columns=physical_units.columns)
            to_output_trace_centered = pd.DataFrame(columns=trace_centered.columns)
            
            for dataset in np.unique(physical_units['dataset']):
                print(dataset)
                for trap_id in np.unique(physical_units[physical_units['dataset'] == dataset]['trap_ID']):
                    for trace in ['A', 'B']:
                        pu_lineage = physical_units[(physical_units['dataset'] == dataset) & (physical_units['trap_ID'] == trap_id) & (physical_units['trace'] == trace)].copy()
                        
                        sampled_lineages = based_on_target(targets, dep, pu_lineage, sg_cov, inter_covariance)
                        
                        to_output_physical_units = to_output_physical_units.append(sampled_lineages, ignore_index=True)
                        
                        # Turn them into trace centered dataframes
                        sampled_lineages[phenotypic_variables] = sampled_lineages[phenotypic_variables] - sampled_lineages[phenotypic_variables].mean()
                        
                        to_output_trace_centered = to_output_trace_centered.append(sampled_lineages, ignore_index=True)
            
            to_output_physical_units.to_csv('{}/{}/{} pu.csv'.format(args.save_folder, args.models_save, filename), index=False)
            to_output_trace_centered.to_csv('{}/{}/{} tc.csv'.format(args.save_folder, args.models_save, filename), index=False)
        exit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process Lineage Data.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                        required=False, default='trace_centered.csv')
    parser.add_argument('-inter_cov', '--inter_cov', metavar='', type=str, help='What to name the intergenerational covariance table.',
                        required=False, default='inter_covariance.csv')
    parser.add_argument('-models_save', '--models_save', metavar='', type=str, help='What to name the intergenerational covariance table.',
                        required=False, default='model_lineage_dataframes')
    
    args = parser.parse_args()
    
    create_folder('{}/{}'.format(args.save_folder, args.models_save))

    first_time = time.time()
    start_time = first_time
    
    main(args)

    print("--- %s seconds ---" % (time.time() - start_time))
