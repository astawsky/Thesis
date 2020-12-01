#!/usr/bin/env bash

import os
import argparse
from CustomFuncsAndVars.global_variables import create_folder
from random import seed
from DataframeCreation.lineage_data_processing import main as lineage_data
from DataframeCreation.pair_data_processing import main as pair_data
from DataframeCreation.intergenerational_correlations import main as inter
from DataframeCreation.pair_gen_specific import main as pair_gen_specific
from DataframeCreation.figure1_analysis import main as figure1_analysis
from DataframeCreation.figure2_analysis import main as figure2_analysis
from DataframeCreation.figure3_analysis import main as figure3_analysis
import time


def main():

    start_time = time.time()
    # --sp_infiles, --nc_infiles, --save_folder, --pu, --tc
    
    lineage_data(args)
    
    print('past lineage data')
    
    # --puc, --tac, --tcc, --number_of_control_lineage_pairs, --random_seed
    
    # So we can compare the trace-centered and the physical units
    seed(args.random_seed)
    
    pair_data(args)
    
    print('past pair data')
    
    # We use --bs and --inter
    
    inter(args)
    
    print('inter')
    
    # --gen_specific_amount
    
    pair_gen_specific(args)
    
    print('par gen specific')
    
    # --population_sampled, --ta, --ebp, --kld
    
    figure1_analysis(args)
    
    print('figure 1 analyses. Namely, KL divergences ')
    
    # --shuffled, --cta, --vcta
    
    figure2_analysis(args)
    
    print('figure 2 analyses.')
    
    # --pair_kld, --lin_and_time
    
    figure3_analysis(args)
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
def debug():
    figure2_analysis(args)
    

parser = argparse.ArgumentParser(description='Process Lineage Data.')
parser.add_argument('-SL', '--sp_infiles', metavar='', type=str,
                    help='Location of Sister Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/SL')
parser.add_argument('-NL', '--nc_infiles', metavar='', type=str,
                    help='Location of Neighboring Cell Pair Raw Data', required=False, default=r'/Users/alestawsky/PycharmProjects/Thesis/RawData/NL')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Data')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                    required=False, default='physical_units.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                    required=False, default='trace_centered.csv')
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
parser.add_argument('-bs', '--bs', metavar='', type=int, help='How many bootstraps per covariance should be done?',
                    required=False, default=0)
parser.add_argument('-inter', '--inter', metavar='', type=str, help='What to name the intergenerational correlation dataframe.',
                    required=False, default='total_inter_correlations.csv')
parser.add_argument('-gsa', '--gen_specific_amount', metavar='', type=int, help='How many generations to calculate for the generation-specific correlations.',
                    required=False, default=6)
parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                    required=False, default='population_lineages.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                    required=False, default='time_averages.csv')
parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                    required=False, default='ergodicity_breaking_parameter.csv')
parser.add_argument('-kld', '--kld', metavar='', type=str,
                    help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                    required=False, default='kullback_leibler_divergences.csv')

parser.add_argument('-shuffled', '--shuffled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the generation-shuffled lineages.',
                    required=False, default='shuffled_generations.csv')
parser.add_argument('-cta', '--cta', metavar='', type=str,
                    help='The filename of the dataframe that contains the expanding time-averages of each kind of lineage (Trace, Population, Shuffled or Trace-centered).',
                    required=False, default='expanding_mean.csv')
parser.add_argument('-vcta', '--vcta', metavar='', type=str,
                    help='The filename of the dataframe that contains the variation of the expanding time-averages of each kind of lineage (Trace, Population, Shuffled or Trace-centered).',
                    required=False, default='variation_of_expanding_mean.csv')
parser.add_argument('-cum_sum', '--cum_sum', metavar='', type=str,
                    help='The filename of the dataframe that contains the cumulative sums (random walks) of each kind of lineage (Trace, Population, Shuffled or Trace-centered).',
                    required=False, default='cumulative_sum.csv')
parser.add_argument('-v_cum_sum', '--v_cum_sum', metavar='', type=str,
                    help='The filename of the dataframe that contains the variation of the cumulative sums (random walks) of each kind of lineage (Trace, Population, Shuffled or Trace-centered).',
                    required=False, default='variation_of_cumulative_sum.csv')
parser.add_argument('-pair_kld', '--pair_kld', metavar='', type=str, help='The filename of the dataframe that contains the kullback-leibler divergences between distributions of all pair lineages.',
                    required=False, default='kullback_leibler_divergences_for_pair_lineages.csv')
parser.add_argument('-lat', '--lin_and_time', metavar='', type=str, help='The filename of the dataframe that contains the model between distributions of all pair lineages over lineages and time.',
                    required=False, default='over_lineages_and_time.csv')

args = parser.parse_args()

create_folder(args.save_folder)

debug()

