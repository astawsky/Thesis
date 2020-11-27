#!/usr/bin/env bash

import os
import argparse
import sys
from CustomFuncsAndVars.global_variables import create_folder
from random import seed
from DataframeCreation.lineage_data_processing import main as lineage_data
from DataframeCreation.pair_data_processing import main as pair_data
from DataframeCreation.intergenerational_correlations import main as inter
import time


first_time = time.time()
start_time = first_time

parser = argparse.ArgumentParser(description='Process Lineage Data.')
parser.add_argument('-SP', '--sp_infiles', metavar='', type=str,
                    help='Location of Sister Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/SP')
parser.add_argument('-NC', '--nc_infiles', metavar='', type=str,
                    help='Location of Neighboring Cell Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/NC')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                    required=False, default='physical_units.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                    required=False, default='trace_centered.csv')

args = parser.parse_args()

print(type(args), args)

create_folder(args.save_folder)

# lineage_data(args)

print('past lineage data')


parser.add_argument('-puc', '--puc', metavar='', type=str, help='What to name the physical units dataframe with control added',
                    required=False, default='physical_units_with_control.csv')
parser.add_argument('-tac', '--tac', metavar='', type=str, help='What to name the time-averages dataframe with control added',
                    required=False, default='time_averages_with_control.csv')
parser.add_argument('-tcc', '--tcc', metavar='', type=str, help='What to name the trace-centered dataframe with control added',
                    required=False, default='trace_centered_with_control.csv')
parser.add_argument('-ctrl_number', '--number_of_control_lineage_pairs', metavar='', type=int, help='How many Control lineage pairs should we make.',
                    required=False, default=85)
parser.add_argument('-seed', '--random_seed', metavar='', type=int, help='Seed for python randomness generator.',
                    required=False, default=42)

args = parser.parse_args()

print(type(args), args)

# So we can compare the trace-centered and the physical units
seed(args.random_seed)

# pair_data(args)

print('past pair data')


parser.add_argument('-bs', '--bs', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)
parser.add_argument('-inter', '--inter', metavar='', type=str, help='What to name the intergenerational correlation dataframe.',
                    required=False, default='total_inter_correlations.csv')

args = parser.parse_args()

print(type(args), args)

inter(args)


parser.add_argument('-gsa', '--gen_specific_amount', metavar='', type=int, help='How many generations to calculate for the generation-specific correlations.',
                    required=False, default=6)


print("--- %s seconds ---" % (time.time() - start_time))
