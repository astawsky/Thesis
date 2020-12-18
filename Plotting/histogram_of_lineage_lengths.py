#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import symbols, units, mm_data_names, create_folder
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main(args):
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.pu))
    
    lengths = []
    for lineage_id in physical_units.lineage_ID.unique():
        trace = physical_units[physical_units['lineage_ID'] == lineage_id].copy()
        lengths.append(len(trace))

    # stylistic reasons
    sns.set_context('paper')
    sns.set_style("ticks", {'axes.grid': True})
    
    # fig, ax = plt.subplots(tight_layout=True, figsize=[3, 3])
    
    sns.displot(data=lengths, label=r'${}$ lineages'.format(len(lengths)) + '\n' + r'$\sim {} \pm {}$ long'.format(np.int(np.mean(lengths)), np.int(np.std(lengths))))
    plt.xlabel('lineage lengths')
    plt.ylabel('PDF')
    plt.title(args.data_origin)
    plt.legend(title='')
    plt.tight_layout()
    plt.savefig('{}/Lineage Length Histogram pdf.png'.format(args.figs_location), dpi=300)
    plt.close()

    sns.displot(data=lengths, label=r'${}$ lineages'.format(len(lengths)) + '\n' + r'$\sim {} \pm {}$ long'.format(np.int(np.mean(lengths)), np.int(np.std(lengths))), kind="ecdf")
    plt.xlabel('lineage lengths')
    plt.ylabel('PDF')
    plt.title(args.data_origin)
    plt.legend(title='')
    plt.tight_layout()
    plt.savefig('{}/Lineage Length Histogram cdf.png'.format(args.figs_location), dpi=300)
    plt.close()


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # How long does running this take?
    first_time = time.time()
    
    for data_origin in mm_data_names + ['SM']:
        parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-pu', '--pu', metavar='', type=str, help='What to name the physical units dataframe.',
                            required=False, default='physical_units.csv')
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        main(args)
