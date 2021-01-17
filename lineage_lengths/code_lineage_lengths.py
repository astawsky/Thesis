#!/usr/bin/env bash

from AnalysisCode.global_variables import symbols, units, dataset_names, create_folder, wang_datasets
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def main(args):
    physical_units = pd.read_csv(args['pu'])
    
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
    plt.title(args['data_origin'])
    plt.legend(title='')
    plt.tight_layout()
    plt.savefig('PDF/{}.png'.format(args['data_origin']), dpi=300)
    plt.close()
    
    sns.displot(data=lengths, label=r'${}$ lineages'.format(len(lengths)) + '\n' + r'$\sim {} \pm {}$ long'.format(np.int(np.mean(lengths)), np.int(np.std(lengths))), kind="ecdf")
    plt.xlabel('lineage lengths')
    plt.ylabel('PDF')
    plt.title(args['data_origin'])
    plt.legend(title='')
    plt.tight_layout()
    plt.savefig('CDF/{}.png'.format(args['data_origin']), dpi=300)
    plt.close()


if __name__ == '__main__':
    import argparse
    import os
    
    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')
    
    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    
    # Finalize the arguments
    input_args = parser.parse_args()
    
    create_folder('PDF')
    create_folder('CDF')
    
    for data_origin in dataset_names:  # input_args.dataset_names:
        
        print(data_origin)
        
        processed_data = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/'
        
        """
                data_origin ==> Name of the dataset we are analysing
                raw_data ==> Where the folder containing the raw data for this dataset is
                processed_data ==> The folder we will put the processed data in
                """
        args = {
            'data_origin': data_origin,
            'pu': processed_data + 'z_score_under_3/physical_units_without_outliers.csv' if data_origin in wang_datasets else processed_data + 'physical_units.csv'
        }
        
        main(args)
