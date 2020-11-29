#!/usr/bin/env bash

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import argparse
from itertools import combinations
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols


def main(args):
    
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
    
    sns.barplot(data=df.replace(symbols['physical_units']), hue='kind', x='variable', y='value', ci=5)
    plt.ylabel('')
    plt.xlabel('')
    plt.show()
    plt.close()
    
    
def other(args):
    
    possible_dependencies = ['previous, ' + variable for variable in ['generationtime', 'growth_rate', 'division_ratio']] + ['length_birth']
    
    actual_dependencies = []
    for n in np.arange(1, len(possible_dependencies)):
        actual_dependencies += [list(val) for val in list(combinations(possible_dependencies, n))]
    
    sg_cov = trace_centered[['generationtime', 'growth_rate', 'length_birth', 'division_ratio']].cov()
    for targets in [['generationtime'], ['growth_rate'], ['generationtime', 'growth_rate']]:
        for dep in actual_dependencies[2:]:
            '{}/{}/{}.csv'.format(args.save_folder, args.models_save, filename)
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process Lineage Data.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
    parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                        required=False, default='physical_units.csv')

    args = parser.parse_args()

    main(args)
