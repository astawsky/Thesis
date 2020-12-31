#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress


def main(args):
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu)).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)

    population_sampled = shuffle_info(physical_units, mm=args.MM)

    lineage_shuffled = shuffle_lineage_generations(physical_units, args.MM)
    
    for lin_id in physical_units.lineage_ID.unique():
        stand_lineage = physical_units[physical_units['lineage_ID'] == lin_id].copy().sort_values('generation')
        stand_cumsum_lineage = ((stand_lineage[phenotypic_variables] - stand_lineage[phenotypic_variables].mean()) / stand_lineage[phenotypic_variables].std()).cumsum()
        for variable in phenotypic_variables:
            lin_walk = stand_cumsum_lineage[variable] - stand_lineage[stand_lineage['generation'] == stand_lineage.generation.min()]
            


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # The final_generations for the different datasets
    fg = {
        'SM': 45,
        # 'lambda_LB': 60,
        'Maryam_LongTraces': 45,
        'MG1655_inLB_LongTraces': 200,
        # 'LAC_M9': np.nan
        '4-28-2017': 30
    }
    
    # How long does running this take?
    first_time = time.time()
    
    total_ebp_dict = {}
    
    print(dataset_names)
    
    # Do all the Mother Machine data
    for data_origin in dataset_names:  # dataset_names
        data_origin = 'lambda_LB'
        print(data_origin)
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ProcessedData/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-pop', '--population_sampled', metavar='', type=str, help='The filename of the dataframe that contains the physical units of the population sampled lineages.',
                            required=False, default='artificial_lineages.csv')
        parser.add_argument('-ta', '--ta', metavar='', type=str, help='What to name the time-averages dataframe.',
                            required=False, default='time_averages.csv')
        parser.add_argument('-ebp', '--ebp', metavar='', type=str, help='What to name the dataframe containing the ergodicity breaking parameter for each variable.',
                            required=False, default='gamma_ta_corrs_per_gen.csv')
        parser.add_argument('-kld', '--kld', metavar='', type=str,
                            help='What to name the dataframe containing the kullback leibler diverges for each variable between the population ensemble and physical units of lineages.',
                            required=False, default='kullback_leibler_divergences.csv')
        parser.add_argument('-MM', '--MM', metavar='', type=bool, help='Is this MM data?', required=False, default=True)
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        folder_name = args.figs_location + '/plot_ebp_per_gen'
        create_folder(args.figs_location)
        
        # folder_name = args.os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/plot_ebp_per_gen_all_datasets'
        
        create_folder(folder_name)
        
        # if data_origin == 'SM':
        #     main(args, first_generation=1, final_generation=25)  # fg[data_origin]
        # else:
        #     main(args, first_generation=0, final_generation=24)
        
        # ebp = pd.read_csv('{}/{}'.format(args.save_folder, args.ebp))
        # total_length_ebp = pd.read_csv('gamma_ta_corrs_{}.csv'.format(args.data_origin))
        # ebp_pergen = pd.read_csv('gamma_ta_corrs_per_gen_{}.csv'.format(args.data_origin))
        
        # # total_ebp_dict.update({args.data_origin: ebp_pergen})
        #
        # # print('moving on to plotting')
        # # plot_ebp_per_gen({args.data_origin: ebp_pergen}, folder_name)
        #
        # # print(total_length_ebp)
        
        total = pd.read_csv('gamma_ta_corrs_{}.csv'.format(args.data_origin))
        
        plot_ebp(total[total['label'] != 'Shuffled'], symbols['time_averages'])
        
        # # set a style on seaborn module
        # sns.set_context('paper')
        # sns.set_style("ticks", {'axes.grid': True})
        # _, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
        # plot_ebp(, ax, symbols['time_averages'])
        # # save the figure
        # plt.savefig('{}/EBP.png'.format(args.figs_location), dpi=300)
        # # plt.show()
        # plt.close()
        
        print('*' * 200)
