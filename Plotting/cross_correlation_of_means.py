#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, dataset_names, symbols, seaborn_preamble, shuffle_info, create_folder, cmap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, pearsonr


def main(args):
    # the parameters of the scaling analysis
    min_ws = 3
    max_ws = None
    window_size_steps = 1

    # import the labeled measured bacteria in physical units
    physical_units = pd.read_csv('{}/physical_units_with_outliers.csv'.format(args.save_folder))
    artificial_lineages = shuffle_info(physical_units, mm=True)

    # Where we will keep the Normalized Cross Correlations
    ncc = pd.DataFrame(columns=['lineage_ID', 'dataset', 'kind', 'variable1', 'variable2', 'lag', 'r'])

    for lin_id in physical_units.lineage_ID.unique():
        print(lin_id)

        # The lineages that come from different "datasets"
        trace = physical_units[physical_units['lineage_ID'] == lin_id].dropna().sort_values('generation')
        shuffled = trace.sample(frac=1).copy()
        artificial = artificial_lineages[artificial_lineages['lineage_ID'] == lin_id].dropna().sort_values('generation')  # .reset_indices(drop=True)
        # white_noise = pd.Series(np.random.normal(0, 1, len(trace)))
        # noisy_wave = pd.Series([np.sin()])

        if len(trace) < 200:
            print('len(trace) < 120')
            continue

        # This is chosen so that we have at least 10 points to get the pearson correlation of the lagged values
        max_lag = min(len(trace) - 10, 100)

        print('max lag', max_lag)

        # This is to speed it up and so we get clear results from the beginning and we sample less after 50 generations
        if max_lag < 50:
            lags = np.arange(0, max_lag, 1)
        else:
            lags = np.arange(0, 50, 1)
            lags = np.append(lags, np.arange(51, max_lag, 3))

        lags = np.arange(0, 200, 8)

        # If the lineage is long enough to have five lags in it at least
        if len(trace) < 15:  # len(trace) < 3 * min_ws:
            print('lineage not long enough')
            continue

        # What type of lineages do we want to look at?
        names_types_of_lineages = ["Trace", "Artificial", "Shuffled"]  # , "White Noise"

        for variable1 in phenotypic_variables:
            print('variable1:', variable1)
            for variable2 in phenotypic_variables:
                # print('variable2', variable2)

                for dataset in names_types_of_lineages:
                    # print('dataset', dataset)

                    # Decide what type of lineage we will look at
                    if dataset == 'Trace':
                        lineage1 = trace[variable1].values
                        lineage2 = trace[variable2].values
                    elif dataset == 'Artificial':
                        lineage1 = artificial[variable1].values
                        lineage2 = artificial[variable2].values
                    elif dataset == 'Shuffled':
                        lineage1 = shuffled[variable1].values
                        lineage2 = shuffled[variable2].values
                    else:
                        # The white noise then
                        lineage1 = np.random.normal(0, 1, len(trace))
                        lineage2 = np.random.normal(0, 1, len(trace))

                    # We will treat both directions of the lag as different results
                    for lag in lags:

                        # This is because [:-0] does not work for np arrays
                        if lag == 0:
                            x, y = lineage1, lineage2
                        else:
                            x, y = lineage1[lag:], lineage2[:-lag]

                        # # Get the normalized cross correlation
                        # r = pearsonr(x, y)[0]
                        #
                        # # Add it to the dataframe
                        # ncc = ncc.append({
                        #     'lineage_ID': lin_id, 'dataset': dataset, 'kind': 'pu_values', 'variable1': variable1, 'variable2': variable2, 'lag': lag, 'r': r
                        # }, ignore_index=True)

                        # We take the cross correlation of the normalized cumulative sum that is supposed to give us the trend
                        lineage_trend1 = np.cumsum((x - np.mean(x)) / np.std(x))
                        lineage_trend2 = np.cumsum((y - np.mean(y)) / np.std(y))

                        # Get the normalized cross correlation
                        r = pearsonr(lineage_trend1, lineage_trend2)[0]

                        # Add it to the dataframe
                        ncc = ncc.append({
                            'lineage_ID': lin_id, 'dataset': dataset, 'kind': 'normalized_cumsum', 'variable1': variable1, 'variable2': variable2, 'lag': lag, 'r': r
                        }, ignore_index=True)

                    # It can't have any NaNs
                    assert not ncc.isnull().values.any()

    # Save them
    ncc.to_csv('{}/normalized_auto_and_cross_correlations_cumsum.csv'.format(args.save_folder), index=False)

    print('finished saving them')
    
    
def plot(args):
    
    # Create the folder we will save
    folder_fig = '{}/normalized_cross_correlation_cumsum'.format(args.figs_location)
    create_folder(folder_fig)
    
    ncc = pd.read_csv('{}/normalized_auto_and_cross_correlations_cumsum.csv'.format(args.save_folder))
    
    seaborn_preamble()
    
    for variable1 in phenotypic_variables:
        print(variable1)
        for variable2 in phenotypic_variables:
            for kind in ncc.kind.unique():
                # for dataset in ncc.dataset.unique():
                #     print(dataset)
                #     if dataset == 'Trace':
                #         color = 'blue'
                #     else:
                #         color = 'orange'
                #
                #     print(variable2)
                #
                #     relevant = ncc[(ncc['variable1'] == variable1) & (ncc['variable2'] == variable2) & (ncc['dataset'] == dataset) & (ncc['kind'] == kind)].copy()
                #
                #     for lin_id in relevant.lineage_ID.unique():
                #         what_to = relevant[(relevant['lineage_ID'] == lin_id)]
                #         plt.plot(what_to.lag.values, what_to.r.values, marker='x', color=color, alpha=1 if kind == 'pu_values' else 0.5)
                
                fig, ax = plt.subplots(tight_layout=True, figsize=[10, 7])

                relevant1 = ncc[(ncc['variable1'] == variable1) & (ncc['variable2'] == variable2) & (ncc['kind'] == kind)].copy()
                
                sns.boxplot(data=relevant1, hue='dataset', x='lag', y='r', showfliers=False)
                
                # # regress the all the lineages
                # slope, intercept, _, _, std_err = linregress(np.log(relevant.lag.values + 1), np.log(relevant.r.values))
                #
                # domain = relevant.lag.unique()
                #
                # plt.plot(domain, [np.exp(intercept) * (d ** slope) for d in domain], label=dataset, ls='--', color=cmap[0 if dataset == 'Trace' else 1])
                plt.axhline(0, color='black')
        
                plt.ylabel(r'$\rho_n(${}$, ${}$)$'.format(symbols['physical_units'][variable1], symbols['physical_units'][variable2]))
                plt.legend(title='')
                # plt.xlim([0, 20])
                plt.ylim([-1, 1])
                plt.savefig('{}/{} {}.png'.format(folder_fig, variable1, variable2), dpi=300)
                # plt.show()
                plt.close()
    
    
def plot_individuals(args):
    
    # Create the folder we will save
    folder_fig = '{}/normalized_cross_correlation_cumsum'.format(args.figs_location)
    create_folder(folder_fig)
    
    ncc = pd.read_csv('{}/normalized_auto_and_cross_correlations_cumsum.csv'.format(args.save_folder))
    
    seaborn_preamble()
    
    for variable1 in phenotypic_variables:
        print(variable1)
        for variable2 in phenotypic_variables:
            for kind in ncc.kind.unique():

                fig, ax = plt.subplots(tight_layout=True, figsize=[10, 7])
                
                for dataset in ncc.dataset.unique():
                    print(dataset)
                    if dataset == 'Trace':
                        color = cmap[0]
                    elif dataset == 'Artificial':
                        color = cmap[1]
                    else:
                        color = cmap[2]

                    print(variable2)

                    relevant = ncc[(ncc['variable1'] == variable1) & (ncc['variable2'] == variable2) & (ncc['dataset'] == dataset) & (ncc['kind'] == kind)].sort_values('lag').copy()

                    for lin_id in relevant.lineage_ID.unique()[:4]:
                        what_to = relevant[(relevant['lineage_ID'] == lin_id)]
                        plt.plot(what_to.lag.values, what_to.r.values, marker='x', color=color, alpha=1)  # 1 if kind == 'pu_values' else 0.5
                        
                plt.axhline(0, color='black')
                plt.ylabel(r'$\rho_n(${}$, ${}$)$'.format(symbols['physical_units'][variable1], symbols['physical_units'][variable2]))
                plt.ylim([-1, 1])
                plt.legend(title='')
                # plt.xlim([0, 20])
                # plt.savefig('{}/{} {}.png'.format(folder_fig, variable1, variable2), dpi=300)
                plt.show()
                plt.close()


if __name__ == '__main__':
    import argparse
    import os
    import time
    
    # pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/MG1655_inLB_LongTraces/physical_units_with_outliers.csv')
    # for lin_id in pu.lineage_ID.unique():
    #     fig, axes = plt.subplots(4, 1, tight_layout=True)
    #     lineage = pu[pu['lineage_ID'] == lin_id]
    #     lineage = (lineage - lineage.mean()) / lineage.std()
    #     axes[0].axhline(0, color='black')
    #     axes[1].axhline(0, color='black')
    #     axes[2].axhline(0, color='black')
    #     axes[3].axhline(0, color='black')
    #     axes[0].plot(lineage['generation'], lineage['fold_growth'])
    #     axes[1].plot(lineage['generation'], lineage['generationtime'])
    #     axes[2].plot(lineage['generation'], lineage['growth_rate'])
    #     axes[3].plot(lineage['generation'], lineage['length_birth'])
    #     axes[0].set_ylabel(symbols['physical_units']['fold_growth'])
    #     axes[1].set_ylabel(symbols['physical_units']['generationtime'])
    #     axes[2].set_ylabel(symbols['physical_units']['growth_rate'])
    #     axes[3].set_ylabel(symbols['physical_units']['length_birth'])
    #     plt.xlabel('generation')
    #     plt.show()
    #     plt.close()
    #
    # exit()
    
    # How long does running this take?
    first_time = time.time()
    
    # Do all the Mother Machine data
    for data_origin in dataset_names:
        data_origin = 'MG1655_inLB_LongTraces'
        print(data_origin)
        
        # if data_origin == 'lambda_LB' or data_origin == 'MG1655_inLB_LongTraces':
        #     continue
        
        parser = argparse.ArgumentParser(description='Process Mother Machine Lineage Data.')
        parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
        parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
        parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                            required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
        args = parser.parse_args()
        
        # main(args)

        plot(args)

        # plot_individuals(args)

        # df = pd.read_csv('{}/scaling_of_ebp_and_cv.csv'.format(args.save_folder))
        #
        # plot_histograms_of_scaling_exponents(args.figs_location, df)
        
        print('*' * 200)
        exit()
    
    # How long did it take to do the mother machine?
    mm_time = time.time() - first_time
    
    data_origin = 'SM'
    
    print(data_origin)
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    args = parser.parse_args()
    
    main(args)
    
    df = pd.read_csv('{}/scaling_of_ebp_and_cv.csv'.format(args.save_folder))
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
