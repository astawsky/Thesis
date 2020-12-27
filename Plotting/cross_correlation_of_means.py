#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, mm_data_names, symbols, seaborn_preamble, shuffle_info, create_folder, cmap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, pearsonr


def plot_histograms_of_scaling_exponents(figure_folder, histogram_of_regression):
    seaborn_preamble()
    
    # Histogram of H of every kind
    for kind in histogram_of_regression.kind.unique():
        sns.boxplot(data=histogram_of_regression[histogram_of_regression['kind'] == kind], x='variable', y='slope', hue='dataset', order=symbols['physical_units'].values(), showfliers=False)
        plt.legend(title='')
        plt.title(kind)
        plt.ylabel('exponent')
        plt.xlabel('')
        plt.savefig('{}/{} per variable.png'.format(figure_folder, kind), dpi=300)
        # plt.show()
        plt.close()


def recreate_loglog_plots(loglog_plot_recreation, figure_folder, individuals=False):
    # Create the folders where we will save the figures
    if individuals:
        create_folder('{}/Individual Scalings'.format(figure_folder))
        for kind in loglog_plot_recreation.kind.unique():
            create_folder('{}/Individual Scalings/{}'.format(figure_folder, kind))
            for dataset in loglog_plot_recreation.dataset.unique():
                create_folder('{}/Individual Scalings/{}/{}'.format(figure_folder, kind, dataset))
    else:
        create_folder('{}/Together Scalings'.format(figure_folder))
        for kind in loglog_plot_recreation.kind.unique():
            create_folder('{}/Together Scalings/{}'.format(figure_folder, kind))
    
    # For every kind of loglog plot
    for kind in loglog_plot_recreation.kind.unique():
        print(kind)
        for variable in loglog_plot_recreation.variable.unique():
            if individuals:
                for dataset in loglog_plot_recreation.dataset.unique():
                    for lin_id in loglog_plot_recreation.lineage_ID.unique():
                        relevant = loglog_plot_recreation[
                            (loglog_plot_recreation['kind'] == kind) & (loglog_plot_recreation['lineage_ID'] == lin_id) & (loglog_plot_recreation['variable'] == variable) & (
                                    loglog_plot_recreation['dataset'] == dataset)].copy()
                        
                        fig, ax = plt.subplots(tight_layout=True)
                        plt.scatter(relevant['window_sizes'].values, relevant['y_to_regress'].values)
                        plt.plot(relevant['window_sizes'].values, relevant['best_fit_line'].values,
                                 label=r'log(y)$= {:.2} + {:.2} \,$log(x)'.format(relevant.intercept.unique()[0], relevant.slope.unique()[0]))
                        plt.ylabel('Variation')
                        plt.xlabel('Window Size')
                        plt.yscale('log')
                        plt.xscale('log')
                        plt.legend(title='')
                        plt.tight_layout()
                        plt.show()
                        # plt.savefig('{}/Individual Scalings/{}/{}/{}.png'.format(figure_folder, kind, dataset, lin_id), dpi=300)
                        plt.close()
                        exit()
            else:
                relevant = loglog_plot_recreation[
                    (loglog_plot_recreation['kind'] == kind) & (loglog_plot_recreation['variable'] == variable) & (loglog_plot_recreation['dataset'].isin(['Trace', 'Artificial']))].copy()
                
                seaborn_preamble()
                
                set_trace_legend, set_artificial_legend = False, False
                
                wss_array = {'Trace': np.array([]), 'Artificial': np.array([])}
                ytr_array = wss_array.copy()
                
                # For each lineage separately
                for ind in relevant.index:
                    # Convert the window sizes and mean squared error per windows size from the dataframe from a string to an array of integers or floats
                    wss = relevant['window_sizes'].loc[ind].strip('][').split(' ')
                    wss = [int(r.split('\n')[0]) for r in wss if r != '']
                    ytr = relevant['y_to_regress'].loc[ind].strip('][').split(', ')
                    ytr = [float(r) for r in ytr if r != '']
                    
                    # If it is a trace lineage plot it and include it in the legend
                    if relevant.loc[ind]['dataset'] == 'Trace':
                        # The constant artificial color
                        color = cmap[0]
                        
                        # Add the analysis curve to the array of its dataset for the population regression
                        wss_array['Trace'] = np.append(wss_array['Trace'], wss)
                        ytr_array['Trace'] = np.append(ytr_array['Trace'], ytr)
                        
                        # Include it in the legend but not more than once
                        if set_trace_legend:
                            plt.plot(wss, ytr, color=color, marker='x')
                        else:
                            plt.plot(wss, ytr, color=color, marker='x', label='Trace')
                            set_trace_legend = True
                    # If it is a artificial lineage plot it and include it in the legend
                    elif relevant.loc[ind]['dataset'] == 'Artificial':
                        # The constant artificial color
                        color = cmap[1]
                        
                        # Add the analysis curve to the array of its dataset for the population regression
                        wss_array['Artificial'] = np.append(wss_array['Artificial'], wss)
                        ytr_array['Artificial'] = np.append(ytr_array['Artificial'], ytr)
                        
                        # Include it in the legend but not more than once
                        if set_artificial_legend:
                            plt.plot(wss, ytr, color=color, marker='x')
                            pass
                        else:
                            plt.plot(wss, ytr, color=color, marker='x', label='Artificial')
                            set_artificial_legend = True
                    # We do not want to see the white noise
                    else:
                        continue
                
                # Get the linear regression of all the trajectories pooled for each dataset
                slope_trace, intercept_trace, _, _, std_err_trace = linregress(np.log(np.array(wss_array['Trace']).flatten()), np.log(np.array(ytr_array['Trace']).flatten()))
                slope_art, intercept_art, _, _, std_err_art = linregress(np.log(np.array(wss_array['Artificial']).flatten()), np.log(np.array(ytr_array['Artificial']).flatten()))
                
                # Plot the best fit line and its parameters
                plt.plot(np.unique(wss_array['Trace']), [np.exp(intercept_trace) * (l ** slope_trace) for l in np.unique(wss_array['Trace'])], ls='--', color='blue', linewidth=3,
                         label=r'$' + str(np.round(intercept_trace, 2)) + r'n^{' + str(np.round(slope_trace, 2)) + r'\pm' + str(np.round(std_err_trace, 3)) + r'}$')
                plt.plot(np.unique(wss_array['Artificial']), [np.exp(intercept_art) * (l ** slope_art) for l in np.unique(wss_array['Artificial'])], ls='--', color='red', linewidth=3,
                         label=r'$' + str(np.round(intercept_art, 2)) + r'n^{' + str(np.round(slope_art, 2)) + r'\pm' + str(np.round(std_err_art, 3)) + r'}$')
                
                plt.title(variable)
                plt.ylabel(r'$F(n)$')
                plt.xlabel('n')
                plt.yscale('log')
                plt.xscale('log')
                plt.legend(title='')
                plt.tight_layout()
                # plt.show()
                plt.savefig('{}/Together Scalings/{}/{}.png'.format(figure_folder, kind, variable), dpi=300)
                plt.close()
    # plt.scatter(window_sizes, y_to_regress, color='blue')


def ergodicity_breaking_in_lineage(lin_id, dataset, variable, lineage, min_ws, max_ws=None, window_size_steps=5):
    # How many cycles/generations in the lineages
    total_length = len(lineage)
    
    # Maximum window size possible
    if not max_ws:
        max_ws = np.floor(total_length / 3)
    
    # What are the window sizes? In order for it not to take a long time we might skip some
    window_sizes = np.arange(min_ws, max_ws, window_size_steps, dtype=int)
    
    # Where we will keep the dataframe
    save_df = pd.DataFrame(columns=['lineage_ID', 'dataset', 'variable', 'ws', 'length', 'num_of_ind_means', 'ebp', 'cv_means', 'var_means'])
    
    for ws in window_sizes:
        # the indices of the lineage with window size ws and windows one next to the other
        starting_points = np.arange(0, total_length, ws, dtype=int)[:-1]
        
        # The windowed/partial time-series
        partial_time_series = np.array([lineage.iloc[sp:sp + ws].values for sp in starting_points])
        
        # The windowed/partial time-series mean
        ind_means = np.array([np.mean(pts) for pts in partial_time_series])
        means_partial_ts = np.array([[np.mean(pts) for _ in np.arange(ws)] for pts in partial_time_series]).flatten()
        
        # Make it so that it is all one array instead of the windows seperately, which we needed for the mean of the window
        partial_time_series = partial_time_series.flatten()
        
        # Necessary
        assert len(means_partial_ts) == len(partial_time_series)
        assert np.abs(np.mean(ind_means) - np.mean(partial_time_series)) < .000001  # mean is a linear operator
        
        # What we want to calculate and save into the dataframe
        ebp = np.var(means_partial_ts) / np.var(partial_time_series)
        cv_means = np.var(ind_means) / np.abs(np.mean(partial_time_series))
        
        # Add it to the dataframe
        save_df = save_df.append({
            'lineage_ID': lin_id, 'dataset': dataset, 'variable': variable, 'ws': ws, 'length': len(partial_time_series), 'num_of_ind_means': len(ind_means), 'ebp': ebp, 'cv_means': cv_means,
            'var_means': np.var(ind_means)
        }, ignore_index=True)
    
    # Good practice
    save_df = save_df.sort_values('ws').reset_index(drop=True)
    
    # So that we can regress successfully
    assert len(window_sizes) == len(save_df)
    
    x = [float(l) for l in save_df['ws']]
    y = save_df['ebp'].values
    
    # Get the linear regression parameters
    slope, intercept, _, _, std_err = linregress(np.log(x), np.log(y))
    
    # if slope == np.nan:
    #     print('slope is NaN')
    #     print(x)
    #     print(y)
    #     exit()
    
    ebp_df = {
        'lineage_ID': lin_id, 'dataset': dataset, 'variable': variable, 'slope': slope, 'intercept': intercept, 'std_err': std_err, 'kind': 'ebp'
    }
    
    # Get the linear regression parameters
    slope_cv, intercept_cv, _, _, std_err_cv = linregress(np.log(x), np.log(save_df['cv_means'].values))
    
    # if slope_cv == np.nan:
    #     print('slope cv is NaN')
    #     print(x)
    #     print(save_df['cv_means'].values)
    #     exit()
    
    cv_df = {
        'lineage_ID': lin_id, 'dataset': dataset, 'variable': variable, 'slope': slope_cv, 'intercept': intercept_cv, 'std_err': std_err_cv, 'kind': 'cv'
    }
    
    return [save_df, ebp_df, cv_df]


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
        white_noise = pd.Series(np.random.normal(0, 1, len(trace)))
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
                plt.savefig('{}/{} {}.png'.format(folder_fig, variable1, variable2), dpi=300)
                # plt.show()
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
    for data_origin in mm_data_names:
        print(data_origin)
        data_origin = 'MG1655_inLB_LongTraces'
        
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
    
    plot_histograms_of_scaling_exponents(args.figs_location, df)
    
    # How long did it take to do the mother machine?
    sm_time = time.time() - (mm_time + first_time)
    
    print("--- took {} mins in total: {} mins for the MM data and {} mins for the SM data ---".format((time.time() - first_time) / 60, mm_time / 60, sm_time / 60))
