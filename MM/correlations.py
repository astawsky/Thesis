#!/usr/bin/env bash

import pandas as pd
import numpy as np
import scipy.io as io
from scipy.stats import linregress, pearsonr
import glob
import os
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder


def main():
    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/MMG1655_physical_units.csv')

    trace_centered = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/MMG1655_trace_centered.csv')

    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/MMG1655_physical_units.csv')
    
    # print(physical_units[['generationtime', 'growth_rate']])
    # print(physical_units['division_ratio'])
    # exit()
    
    # for lineage in [physical_units[physical_units['trap_ID'] == trap_id].sort_values('generation') for trap_id in physical_units['trap_ID'].unique()]:
    #     plt.plot(np.log(lineage['length_birth'].values), marker='.', label='measured')
    #     new_length_birth = [lineage[lineage['generation'] == 0]['length_birth'].values[0]]
    #     for gr, itd, div_rat, fold_growth in zip(lineage['growth_rate'].values, lineage['generationtime'].values, lineage[lineage['generation'] > 0]['division_ratio'].values, lineage['fold_growth'].values):
    #         # if div_rat == np.nan:
    #         #     div_rat = .5
    #         new_length_birth.append(new_length_birth[-1] * np.exp(gr*itd) * div_rat)
    #     print(new_length_birth)
    #     plt.plot(np.log(new_length_birth), marker='.', ls='--', label='shuffled')
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.show()
    #     plt.close()
    
    # print(physical_units['length_birth'])
    
    # print(physical_units.corr())
    # sns.heatmap(trace_centered[phenotypic_variables].corr(), center=0, annot=True)
    # plt.show()
    # plt.close()
    
    repeats = []
    for var1 in phenotypic_variables:
        for var2 in phenotypic_variables:
            if (var2 != var1) & (var2 not in repeats):
                indices = physical_units.sort_values(['trap_ID', 'generation']).copy()[[var1, var2]]
                print(indices)
                indices = indices.dropna()
                print(indices)
                print('------')
                # exit()
                
                x = indices[var1].copy()
                y = indices[var2].copy()
                print(x, y)
                print(linregress(x, y))
                slope, intercept = linregress(x, y)[:2]
                plt.scatter(x, y)
                plt.plot(x, [intercept + point * slope for point in x], label='pop: {:.2}'.format(slope))
                plt.xlabel(var1)
                plt.ylabel(var2)
                plt.legend()
                plt.show()
                plt.close()
        repeats.append(var1)


if __name__ == '__main__':
    main()
