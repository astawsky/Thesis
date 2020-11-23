#!/usr/bin/env bash

import sys
import os
import pandas as pd
import numpy as np
import argparse

sys.path.append(r'/Users/alestawsky/PycharmProjects/Thesis')
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder, datasets


parser = argparse.ArgumentParser(description='Dataframes containing: KL divergences, Population physical_units lineages, and the ergodicity breaking parameter for both kinds of lineages.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/physical_units_with_control.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/time_averages_with_control.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/trace_centered_with_control.csv')
parser.add_argument('-bs', '--bootstraps', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)

args = parser.parse_args()

create_folder(args.save_folder)

# save it to the Data folder
pair_correlation = pd.read_csv('{}/over_lineages.csv'.format(args.save_folder))

print(pair_correlation)
print(pair_correlation.columns)

for var_a in phenotypic_variables:
    for var_b in phenotypic_variables:

