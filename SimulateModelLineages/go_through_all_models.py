#!/usr/bin/env bash

# import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
# import glob
# from sklearn.linear_model import LinearRegression
from SimulateModelLineages.create_lineage_dataframes import main as simulate
# from CustomFuncsAndVars.global_variables import phenotypic_variables, datasets, create_folder
from itertools import combinations
import argparse
import os
import time

first_time = time.time()
start_time = first_time

# All possible dependencies
acd_time_array = [[None]]
acd_rate_array = [[None]]
for n in np.arange(1, 5):
    acd_time_array += [list(val) for val in list(combinations(['generationtime', 'growth_rate', 'division_ratio'], n))]
    acd_rate_array += [list(val) for val in list(combinations(['growth_rate', 'division_ratio'], n))]
    
sgd_time_array = [[None]]
sgd_rate_array = [[None], ['generationtime']]
for n in np.arange(1, 3):
    sgd_time_array += [list(val) for val in list(combinations(['growth_rate', 'length_birth'], n))]
    # sgd_rate += [list(val) for val in list(combinations(['generationtime', 'length_birth'], n))]

print(sgd_time_array)
print(sgd_rate_array)
print(acd_time_array)
print(acd_rate_array)
exit()

# print(sgd_rate, len(sgd_rate))
# print(sgd_time, len(sgd_time))
# exit()
# print(acd_time)
# print(len(acd_time))
# print(acd_rate)
# print(len(acd_rate))

parser = argparse.ArgumentParser(description='Process Lineage Data.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                    required=False, default='physical_units.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                    required=False, default='trace_centered.csv')
parser.add_argument('-acd_time', '--acd_time', metavar='', type=str, nargs="+", action='append', help='What auto-correlation dependencies does inter-division time have?',
                    required=False, default=acd_time_array)
parser.add_argument('-acd_rate', '--acd_rate', metavar='', type=str, nargs="+", action='append', help='What auto-correlation dependencies does growth rate have?',
                    required=False, default=acd_rate_array)
parser.add_argument('-sgd_time', '--sgd_time', metavar='', type=str, nargs="+", action='append', help='What same-generation dependencies does inter-division time have?',
                    required=False, default=sgd_time_array)
parser.add_argument('-sgd_rate', '--sgd_rate', metavar='', type=str, nargs="+", action='append', help='What same-generation dependencies does growth rate have?',
                    required=False, default=sgd_rate_array)
parser.add_argument('-inter_cov', '--inter_cov', metavar='', type=str, help='What to name the intergenerational covariance table.',
                    required=False, default='inter_covariance.csv')

args = parser.parse_args()

# What to pass to create_lineage_dataframes.py
auto_corr_dependencies = {
    'generationtime': args.acd_time,
    'growth_rate': args.acd_rate
}

same_gen_dependencies = {
    'generationtime': args.sgd_time,
    'growth_rate': args.sgd_rate
}

# # For all combinations models
# for sgd_time in sgd_time_array:
#     for sgd_rate in sgd_rate_array:
#         for acd_time in acd_time_array:
#             for acd_rate in acd_rate_array:
#                 simulate(args, auto_corr_dependencies, same_gen_dependencies)
#
# # Only single-generation mechanisms
# for sgd_time in sgd_time_array:
#     for sgd_rate in sgd_rate_array:
#         pass
#
# # Only auto-correlation mechanisms
# for acd_time in acd_time_array:
#     for acd_rate in acd_rate_array:
#         pass

simulate(args)
