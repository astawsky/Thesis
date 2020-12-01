#!/usr/bin/env bash

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import argparse
from itertools import combinations
import glob
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, seaborn_preamble


def main(args, variables_to_look_at=['length_birth', 'fold_growth', 'division_ratio', 'growth_rate', 'generationtime']):
    symbolim = symbols['old_var'].copy()
    symbolim.update({'same_length_birth': symbols['physical_units']['length_birth']})
    print(symbolim)
    simulation_files = glob.glob(args.simulation_locations + '/*.csv')
    # print(simulation_files)
    # print([sim.split('/')[-1] for sim in simulation_files])
    # print([sim.split('/')[-1].split('.csv')[0] for sim in simulation_files])
    # print([sim.split('/')[-1].split('.csv')[0].split(' and dependent on ') for sim in simulation_files])
    targets = [sim.split('/')[-1].split('.csv')[0].split(' and dependent on ')[0].split('target ')[-1] for sim in simulation_files]
    dependencies = [sim.split('/')[-1].split('.csv')[0].split(' and dependent on ')[1].split(', ') for sim in simulation_files]
    dependencies = [[d[:-3] if d == dep[-1] else d for d in dep] for dep in dependencies]
    print(len(targets))
    print(len(dependencies))
    print(targets)
    print(dependencies)
    
    # df = pd.DataFrame(columns=['variable', 'model', 'std', 'cv', 'dataset', 'trap_ID', 'trace'])
    # df_exp = pd.DataFrame(columns=['variable', 'model', 'expected_std', 'expected_cv'])
    #
    # for sim in simulation_files:
    #     pu_sim = pd.read_csv(sim)
    #     target = sim.split('/')[-1].split('.csv')[0].split(' and dependent on ')[0].split('target ')[-1]
    #     dep = [d[:-3] if d == sim.split('/')[-1].split('.csv')[0].split(' and dependent on ')[1].split(', ')[-1] else d for d in
    #            sim.split('/')[-1].split('.csv')[0].split(' and dependent on ')[1].split(', ')]
    #
    #     model = ''
    #     for trgt in [target]:
    #         if len(trgt.split(', ')) > 1:
    #             for t in trgt.split(', '):
    #                 model += symbolim[t]
    #         else:
    #             model += symbolim[trgt]
    #
    #     model += '|'
    #
    #     for d in dep:
    #         model += symbolim[d]
    #
    #     print(model)
    #
    #     for variable in variables_to_look_at:
    #         print(variable)
    #         for dataset in np.unique(pu_sim['dataset']):
    #             for trap_id in np.unique(pu_sim[pu_sim['dataset'] == dataset]['trap_ID']):
    #                 # print(trap_id)
    #                 for trace in ['A', 'B']:
    #                     pu_lineage = pu_sim[(pu_sim['dataset'] == dataset) & (pu_sim['trap_ID'] == trap_id) & (pu_sim['trace'] == trace)][variables_to_look_at].copy()
    #                     df = df.append(
    #                         pd.DataFrame(
    #                             {
    #                                 'variable': variable,
    #                                 'model': model,
    #                                 'std': pu_lineage[variable].std(),
    #                                 'cv': pu_lineage[variable].std() / pu_lineage[variable].mean(),
    #                                 'dataset': dataset,
    #                                 'trap_ID': trap_id,
    #                                 'trace': trace
    #                             }, index=[len(df)]
    #                         )
    #                     )
    #
    #         interesting = df[(df['variable'] == variable) & (df['model'] == model)].copy()
    #
    #         to_add = {
    #             'variable': variable,
    #             'model': model,
    #             'expected_std': interesting['std'].mean(),
    #             'expected_cv': interesting['cv'].mean()
    #         }
    #
    #         df_exp = df_exp.append(to_add, ignore_index=True)
    #
    # print('\n'*4)
    # print(df)
    # print(df_exp)
    #
    # df.to_csv('{}/models std and cv.csv'.format(args.save_folder), index=False)
    # df_exp.to_csv('{}/models and expected std and cv.csv'.format(args.save_folder), index=False)
    # print('saved both')
    
    df = pd.read_csv('{}/models std and cv.csv'.format(args.save_folder))
    df_exp = pd.read_csv('{}/models and expected std and cv.csv'.format(args.save_folder))
    
    seaborn_preamble()
    # sns.boxplot(data=df, x='model', y='cv')  # or 'std' as y
    sns.barplot(data=df_exp, x='model', y='expected_cv')
    plt.grid(True)
    plt.ylabel('Expected Coefficient of Variation over lineages')
    plt.legend(title='')
    plt.xlabel('')
    plt.tight_layout()
    # save the figure
    plt.savefig('{}/Model Expected CVs.png'.format(args.figs_location), dpi=300)

    # plt.show()
    plt.close()
    
    
def now_with_trace(args):
    
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    # fig, axes = plt.subplots()
    
    df = pd.DataFrame(columns=['variable', 'kind', 'value'])
    
    for dataset in np.unique(physical_units['dataset']):
        for trap_id in np.unique(physical_units[physical_units['dataset'] == dataset]['trap_ID']):
            for trace in ['A', 'B']:
                
                pu_lineage = physical_units[(physical_units['dataset'] == dataset) & (physical_units['trap_ID'] == trap_id) & (physical_units['trace'] == trace)][phenotypic_variables].copy()
                to_add = pu_lineage.var() / pu_lineage.mean()
                
                for variable in phenotypic_variables:
                    df = df.append({'variable': variable, 'kind': 'Trace-Centered', 'value': to_add[variable]}, ignore_index=True)
    
    to_add = physical_units[phenotypic_variables].var() / physical_units[phenotypic_variables].mean()
    
    for variable in phenotypic_variables:
        df = df.append({'variable': variable, 'kind': 'total', 'value': to_add[variable]}, ignore_index=True)
    
    sns.barplot(data=df.replace(symbols['physical_units']), hue='kind', x='variable', y='value', ci=100)
    plt.grid(True)
    plt.ylabel('Coefficient of Variation')
    plt.legend(title='')
    plt.xlabel('')
    plt.tight_layout()
    # save the figure
    plt.savefig('{}/Trace vs pop CV.png'.format(args.figs_location), dpi=300)
    
    # plt.show()
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process Lineage Data.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
    parser.add_argument('-sim_loc', '--simulation_locations', metavar='', type=str, help='Where the dataframes of all the simulations are located.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/model_lineage_dataframes')
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures')
    
    args = parser.parse_args()
    
    main(args, variables_to_look_at=['length_birth'])
    # now_with_trace(args)
