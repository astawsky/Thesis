#!/usr/bin/env bash

import os
import argparse
from CustomFuncsAndVars.global_variables import create_folder
from Plotting.intergenerational_figures import main as inter_fig
from Plotting.pair_gen_specific_figures import main as pair_gen_spec
from Plotting.single_cell_correlations import main as single_cell_correlations
from Plotting.growth_mechanism import main as growth_mechanism
from Plotting.figure1 import main as figure1
from Plotting.figure2 import main as figure2
from Plotting.figure3 import main as figure3
import time

first_time = time.time()
start_time = first_time

parser = argparse.ArgumentParser(description='Process Lineage Data.')
parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                    required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Data')
parser.add_argument('-inter_fig', '--inter_fig', metavar='', type=str, help='What to name the intergenerational correlation figure.',
                    required=False, default='intergenerational_big_figure.png')
parser.add_argument('-inter', '--inter', metavar='', type=str, help='What to name the intergenerational correlation dataframe.',
                    required=False, default='total_inter_correlations.csv')
parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                    required=False, default=os.path.dirname(os.path.abspath(__file__)) + '/Figures')
parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                    required=False, default='physical_units.csv')
parser.add_argument('-tc', '--tc', metavar='', type=str, help='What to name the trace-centered dataframe.',
                    required=False, default='trace_centered.csv')
parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                    required=False, default='time_averages.csv')
parser.add_argument('-scc', '--scc', metavar='', type=str, help='Where the single cell correlation figures are saved.',
                    required=False, default='single_cell_correlations')
parser.add_argument('-scch', '--scch', metavar='', type=str, help='Where the single cell correlation heatmap figures are saved.',
                    required=False, default='single_cell_correlations_heatmaps')
parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                    required=False, default='ergodicity_breaking_parameter.csv')
parser.add_argument('-kld', '--kld', metavar='', type=str,
                    help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                    required=False, default='kullback_leibler_divergences.csv')
parser.add_argument('-pair_kld', '--pair_kld', metavar='', type=str, help='The filename of the dataframe that contains the kullback-leibler divergences between distributions of all pair lineages.',
                    required=False, default='kullback_leibler_divergences_for_pair_lineages.csv')
parser.add_argument('-lat', '--lin_and_time', metavar='', type=str, help='The filename of the dataframe that contains the model between distributions of all pair lineages over lineages and time.',
                    required=False, default='over_lineages_and_time.csv')


args = parser.parse_args()

create_folder(args.save_folder)
create_folder(args.figs_location)

inter_fig(args)

pair_gen_spec(args)

single_cell_correlations(args)

growth_mechanism(args)

figure1(args)

figure2(args)

figure3(args)

print("--- %s seconds ---" % (time.time() - start_time))
