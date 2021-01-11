#!/usr/bin/env bash

from AnalysisCode.global_variables import sm_datasets, wang_datasets, dataset_names, create_folder, shuffle_info, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, symbols
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr


# The function that will do the analysis
def get_deviation_dataframes(variables, msd_df, regs, dataframe, label, time_series):
    for count, lin_id in enumerate(dataframe.lineage_ID.unique()):
        # print(count)
        
        # Temporary dictionary to make the process of saving data faster
        temp_dict = {}
        ind = 0
        
        # Get the lineage from the dataframe
        lineage = dataframe[dataframe['lineage_ID'] == lin_id].copy().sort_values('generation')
        
        if len(lineage) < 20:
            # print('{} < 20 => too short for our analysis'.format(len(lineage)))
            continue
            
        distances = np.arange(1, len(lineage) - 15)
        
        # Do we look at the cumulative sum or the signal itself?
        if time_series == 'Signal':
            # Get the cumulative sum for the variables in the lineage
            stand_cumsum_lineage = ((lineage[variables] - lineage[variables].mean()) / lineage[variables].std())
        else:
            # Get the cumulative sum for the variables in the lineage
            stand_cumsum_lineage = ((lineage[variables] - lineage[variables].mean()) / lineage[variables].std()).cumsum()
        
        for variable in variables:
            
            # the random walk
            lin_walk = stand_cumsum_lineage[variable].values - stand_cumsum_lineage[variable].dropna().iloc[0]
            
            msd_lineage = []
            
            for distance in distances:
                
                msd = np.mean(np.array([
                    # The calculation
                    (first - second) ** 2
                    
                    # the points in the time-series
                    for first, second in zip(lin_walk[:-distance], lin_walk[distance:])
                    
                    # Only include the points that are not NaNs
                    if ~np.any([np.isnan(first), np.isnan(second)])
                ]))
                
                msd_lineage.append(msd)
                
                temp_dict[ind] = {
                    'MSD': msd,
                    'RMSD': np.sqrt(msd),
                    'NRMSD': np.sqrt(msd) / (np.nanmax(lin_walk) - np.nanmin(lin_walk)),
                    'CV_NRMSD': np.sqrt(msd) / np.nanmean(lin_walk),
                    'distance': distance,
                    'dataset': args['data_origin'],
                    'n': len(lin_walk) - distance,
                    'variable': variable,
                    'lineage': lin_id,
                    'kind': label
                }
                
                # make space in the dictionary
                ind += 1
                
            # print(temp_dict)
            # exit()
            # print(np.log(distances))
            # print(np.log(msd_lineage))
            # exit()
            
            # Get the slope of the binned data...
            slope, intercept, _, _, std_err = linregress(np.log(distances), np.log(msd_lineage))
            
            regs = regs.append({
                'intercept': np.exp(intercept),
                'slope': slope,
                'std_error': std_err,
                'variable': variable,
                'lineage': lin_id,
                'kind': label
            }, ignore_index=True)
            
            # # Plot it
            # plt.plot(lin_walk)
            # plt.axhline(0, c='black')
            # plt.show()
            # plt.close()
        # print(temp_dict)
        
        # Append them to the big dataframes
        msd_df = msd_df.append(pd.DataFrame.from_dict(temp_dict, "index"), ignore_index=True)
    
    return [msd_df, regs]


def main(args, variables, time_series, folder_to_save_dfs):
    # Import the physical units
    physical_units = pd.read_csv('{}/physical_units.csv'.format(args['processed_data'])).sort_values(['lineage_ID', 'generation']).reset_index(drop=True)
    
    # Create the population sampled dataframes
    population_sampled = shuffle_info(physical_units, mm=args['MM'])
    
    # Create the shuffled generations inside individual lineages dataframe
    lineage_shuffled = shuffle_lineage_generations(physical_units, args['MM'])
    
    # This is where we will store the results of our analysis
    mean_square_deviations = pd.DataFrame(columns=['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD', 'distance', 'dataset', 'n', 'variable', 'lineage', 'kind'])
    
    regressions = pd.DataFrame(columns=['intercept', 'slope', 'std_error', 'variable', 'lineage', 'kind'])
    
    # For each type of lineage, do the analysis
    for dataframe, label in zip([physical_units, population_sampled, lineage_shuffled], ['Trace', 'Artificial', 'Shuffled']):
        print(label)
        
        mean_square_deviations, regressions = get_deviation_dataframes(variables, mean_square_deviations, regressions, dataframe, label, time_series)
    
    # Save them to the folder
    mean_square_deviations.to_csv('{}/MSD_df.csv'.format(folder_to_save_dfs), index=False)
    regressions.to_csv('{}/regression_parameters_df.csv'.format(folder_to_save_dfs), index=False)
    
    # Make sure all the regressions went well
    assert ~regressions.isnull().values.any()
    
    return [mean_square_deviations, regressions]


def plot_individual_msds(args, variables, where_to_save_figures, dataframe, label, loglog, xlim):
    
    # create the folders
    for y in ['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD']:
        if y != 'MSD':
            continue
        # Create the folders where we will save the figures
        create_folder(where_to_save_figures + '/' + y)
        create_folder(where_to_save_figures + '/' + y + '/' + label)
        if xlim != False:
            create_folder(where_to_save_figures + '/' + y + '/' + label + '/cut_at_20')
    
    seaborn_preamble()
    
    for variable in variables:
        # if variable != 'fold_growth':
        #     continue
        
        for y in ['MSD', 'RMSD', 'NRMSD', 'CV_NRMSD']:
            if y != 'MSD':
                continue
            
            relevant = dataframe[(dataframe['variable'] == variable)].sort_values(['lineage', 'distance']).copy()
            
            # slope, intercept, _, _, std_err = linregress(np.log(relevant['distance'].values), np.log(relevant[y].values))
            
            for lin_id in relevant.lineage.unique():
                to_plot = relevant[relevant['lineage'] == lin_id].sort_values('distance')
                plt.plot(to_plot['distance'].values, to_plot[y].values)  # , marker='o' to_plot['distance'].values
            
            # sns.pointplot(data=relevant, x='distance', y='MSD', color='black', ci=95)
            #
            # relevant['distance'] = relevant['distance'] - 1
            # sns.lineplot(data=relevant, x='distance', y=y, ci=None, marker=None, hue='lineage')
            
            the_y_mean = np.array([relevant[relevant['distance'] == dist][y].mean() for dist in relevant.distance.unique()])
            yerr = np.array([relevant[relevant['distance'] == dist][y].std() for dist in relevant.distance.unique()])
            
            plt.errorbar(x=relevant.distance.unique(), y=the_y_mean, yerr=yerr, color='black', marker='o', capsize=1, capthick=1)
            
            # plt.plot(relevant.distance.unique(), [np.exp(intercept) * (distance ** slope) for distance in relevant.distance.unique()], ls='--', linewidth=3,
            #          label=str(np.exp(intercept))[:4] + r'$x^{' + str(slope)[:4] + '}$', color='black')
            
            plt.title(args['data_origin'] + ': ' + label)
            # plt.legend(title='')
            plt.ylabel(r'{}('.format(y) + symbols['physical_units'][variable] + r')')
            plt.xlabel('Distance')
            if xlim != False:
                plt.xlim([1, xlim])
                save_figs = where_to_save_figures + '/' + y + '/' + label + '/cut_at_20'
            else:
                save_figs = where_to_save_figures + '/' + y + '/' + label
            plt.tight_layout()
            
            # log or not
            if loglog:
                plt.yscale('log')
                plt.xscale('log')
                plt.savefig('{}/{}_loglog'.format(save_figs, variable), dpi=300)
            else:
                plt.savefig('{}/{}'.format(save_figs, variable), dpi=300)
                
            # plt.show()
            plt.close()


# def heatmap_analogs(args, df, label, variables=phenotypic_variables[3:], suffix=''):
#     sns.set_context('paper')
#     sns.set_style("ticks", {'axes.grid': True})
#
#     latex_symbols = {variable: symbols[label][variable] for variable in variables}
#
#     df = df[variables].copy().rename(columns=latex_symbols)
#
#     to_plot = pg.pairwise_corr(df, method='pearson')
#
#     # reformat the results into a dataframe we can use for the heatmap
#     to_plot_df = pd.DataFrame(columns=to_plot.X.unique(), index=to_plot.Y.unique(), dtype=float)
#     repeats = []
#     for x in to_plot.X.unique():
#         for y in to_plot.Y.unique():
#             if x != y and y not in repeats:
#                 # Add the results
#                 to_plot_df[x].loc[y] = to_plot[(to_plot['X'] == x) & (to_plot['Y'] == y)]['r'].values[0]
#
#         repeats.append(x)
#
#     # The mask to show only the lower triangle in the heatmap
#     mask = np.ones_like(to_plot_df)
#     mask[np.tril_indices_from(mask)] = False
#
#     fig, ax = plt.subplots(tight_layout=True, figsize=[7 * (2 / 3), 5.5 * (2 / 3)])
#
#     sns.heatmap(to_plot_df, annot=True, square=True, vmin=-1, vmax=1, mask=mask, center=0)
#     # plt.title(label)
#     # plt.tight_layout()
#     # plt.savefig('{}/{}/{}{}.png'.format(args.figs_location, args.scch, label, suffix), dpi=300)
#     plt.show()
#     plt.close()


def hists_of_regression_stats(dataframe, where_to_save_figure):
    name = '{}/Histogram of regression parameters'.format(where_to_save_figure)
    create_folder(name)
    seaborn_preamble()
    df = dataframe.replace(symbols['physical_units']).copy()
    for y in ['slope', 'intercept', 'std_error']:
        sns.boxplot(data=df, x='variable', y=y, showfliers=False, hue='kind')
        plt.title('')
        plt.xlabel('')
        plt.legend(title='')
        plt.tight_layout()
        # plt.show()
        plt.savefig(name + '/' + y + '.png', dpi=300)
        plt.close()


if __name__ == '__main__':
    import argparse
    import os
    import time

    # Create the arguments for this function
    parser = argparse.ArgumentParser(description='Decide which datasets to process Mother Machine and Sister Machine Raw Data for.')

    parser.add_argument('-dataset_names', '--dataset_names', metavar='', nargs="+", help='What is the label for this data for the Data and Figures folders?', required=False,
                        default=dataset_names)
    parser.add_argument('-cumsum', '--cumsum', metavar='', type=bool, help='Do we want to run the MSD analysis for the cumulative sum?', required=False,
                        default=True)
    parser.add_argument('-signal', '--signal', metavar='', type=bool, help='Do we want to run the MSD analysis for the phenotypic variable values in generation-time?', required=False,
                        default=True)

    # Finalize the arguments
    input_args = parser.parse_args()
    
    # How long does running this take?
    first_time = time.time()
    
    types_to_do_analysis = []
    
    if input_args.cumsum:
        types_to_do_analysis.append('Cumsum')
    if input_args.signal:
        types_to_do_analysis.append('Signal')
    if len(types_to_do_analysis) == 0:
        raise IOError('didnt specify type!')
    
    # Do all the Mother Machine data
    for data_origin in input_args.dataset_names:
        print(data_origin)
        if data_origin not in wang_datasets:
            continue

        """
                data_origin ==> Name of the dataset we are analysing
                processed_data ==> The folder we will put the processed data in
                MM ==> If the data comes from the mother machine or the sister machine
        """
        args = {
            'data_origin': data_origin,
            'processed_data': os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Datasets/' + data_origin + '/ProcessedData/',
            'MM': data_origin not in sm_datasets,
            'analysis_dataframes_path': r'/Users/alestawsky/PycharmProjects/Thesis/MSD'
        }

        for time_series in types_to_do_analysis:
            print(time_series)
        
            # create_folder(args.figs_location)
            
            level1 = args['analysis_dataframes_path'] + '/{}'.format(time_series)
            
            create_folder(level1)
            
            level2 = level1 + '/Dataframes'
            
            create_folder(level2)
            
            fig_level2 = level1 + '/Figures'
            
            create_folder(fig_level2)
            
            level3 = level2 + '/' + args['data_origin']
            
            create_folder(level3)
            
            fig_level3 = fig_level2 + '/' + args['data_origin']
            
            create_folder(fig_level3)
            
            if data_origin == 'SM':
                main(args, phenotypic_variables, time_series, level3)
            else:
                main(args, phenotypic_variables, time_series, level3)
                
            print('Finished analysis')
            
            out_df = pd.read_csv('{}/MSD_df.csv'.format(level3))
            reg_df = pd.read_csv('{}/regression_parameters_df.csv'.format(level3))
            
            # Plot the boxplots of the regression parameters
            hists_of_regression_stats(reg_df, fig_level3)
            
            print('Finished Histograms')
            
            # Plot the individual msds for different lineage types and variables
            for label in out_df.kind.unique():
                # print(label)
                
                # only use the trace, artificial or shuffled lineages
                dataframe = out_df[out_df['kind'] == label].copy()
                
                # loglog and normal scale
                for logl in [True, False]:
                    
                    # The ones that are cutoff at 20 and the ones that are not
                    for xlim in [False, 20]:
                        
                        # Plot it
                        plot_individual_msds(args, phenotypic_variables, fig_level3, dataframe, label, logl, xlim)
                
                print('Finished plotting the MSD by generational-distance of '+label+' lineages')
        print('*' * 200)
