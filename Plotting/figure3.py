#!/usr/bin/env bash

from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, datasets, symbols, hierarchy, seaborn_preamble
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
    
    
# def main(args):
#     """ The KL-Divergence """
#
#     def kl_div_plot(ax):
#         # Here we start plotting both KL-Divergences, The good part with this is that it takes into consideration the standard deviation as well
#         # sns.barplot(data=kld_df, x='combo', y='kldiv 1', hue='dataset', hue_order=['S', 'NL', 'C'], order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ci=33,
#         #             ax=ax, edgecolor='black')
#         sns.boxplot(data=kld_df.replace(symbols['physical_units']), x='variable', y='value', hue='dataset', hue_order=datasets, order=hierarchy, showfliers=False, ax=ax, palette=cmap)
#         ax.legend()
#         ax.set_xlabel('')
#         ax.set_ylabel(r'$D_{KL}$')
#         ax.legend()
#
#     """ diagonal pair plots """
#
#     def diagonal_corrs(kind, ax):
#
#         # print(over_time_and_lineages.columns)
#         #
#         # print(np.unique(over_time_and_lineages.AB_parameter_pair))
#         #
#         # # where we store what to plot
#         # plot_df = pd.DataFrame(columns=over_time_and_lineages.columns)
#         #
#         # # # put in the diagonal parameter pair correlations
#         # # for label in np.unique(over_time_and_lineages.AB_parameter_pair):
#         # #     plot_df = plot_df.append(over_time_and_lineages[(over_time_and_lineages['AB_parameter_pair'] == label)], ignore_index=True)
#         #
#         # plot_df = plot_df[plot_df['kind'] == kind].copy()
#
#         sns.barplot(data=over_time_and_lineages[over_time_and_lineages['kind'] == kind], x='AB_parameter_pair', y='correlation', hue='dataset', hue_order=datasets,
#                     order=symbols[kind].values(), ci=95, dodge=True, ax=ax, edgecolor='black', palette=cmap)
#         ax.set_xlabel('')
#         ax.set_ylabel('')
#         ax.get_legend().remove()
#         ax.set_ylim([-.05, .38])
#         # plt.show()
#         # plt.close()
#
#     """ the pearson correlation coefficient between all same-generation lineage pair-bacteria in each dataset separately in Physical Units """
#
#     def diagonal_units(ax, filename='physical_units', n_boots=500):
#         # This is where we will save the dataframe keeping the correlations
#         directory = '{}/{} correlation dataframes, {} bootstraps'.format(filename, n_boots)
#
#         # the one we will plot
#         corr_df = pd.DataFrame(columns=['AB_parameter_pair', 'correlation', 'dataset'])
#
#         # Gather the data into a bigger dataframe
#         for variable in phenotypic_variables:
#             # We will add this pair parameter selection
#             to_add = pd.read_csv(directory + '/{}, {}.csv'.format(variable, variable))
#             print('to_add')
#             print(to_add)
#             to_add['AB_parameter_pair'] = r'{}'.format(symbols['physical_units'][variable], symbols['physical_units'][variable])
#
#             # contains all the correlations of all pair phenotypic_variables
#             corr_df = corr_df.append(to_add, ignore_index=True).reset_index(drop=True)
#
#         print('corr df')
#         print(corr_df)
#
#         # # Stylistic purposes
#         # sns.set_context('paper')
#         # sns.set_style("whitegrid", {'axes.grid': True})
#
#         # plot this dataframe
#         # ax.yaxis.grid(which='major')
#         sns.boxplot(
#             x='AB_parameter_pair',
#             y='correlation',
#             hue='dataset',
#             hue_order=datasets,
#             showfliers=False,
#             order=hierarchy,
#             data=corr_df,
#             ax=ax,
#             palette=cmap
#         )
#         ax.axhline(0, ls='-', color='black')
#         ax.set_xlabel('')
#         ax.set_ylabel(r'$\rho(\cdot, \cdot)$')  # , fontsize='x-large'
#         ax.set_xticklabels(ax.get_xticklabels())  # , fontsize='x-large'
#         # ax.legend(title='') # , fontsize='large'
#         ax.get_legend().remove()
#         # plt.show()
#         # plt.close()
#
#     # read the csv file
#     kld_df = pd.read_csv('{}/{}'.format(args.save_folder, args.pair_kld))  # 'Data/kullback_leibler_divergences_for_pair_lineages.csv'
#     over_time_and_lineages = pd.read_csv('{}/{}'.format(args.save_folder, args.lin_and_time))  # 'Data/bootstrapped pearson correlations of time-averages.csv'
#
#     seaborn_preamble()
#
#     # create the figure and construct the layout of the figure
#     fig, axes = plt.subplots(2, 2, tight_layout=True, figsize=[7, 7])
#     axes = axes.flatten()
#
#     cmap = sns.color_palette('tab10')[5:8]
#
#     # plot the figures
#     kl_div_plot(axes[0])
#     diagonal_corrs('physical_units', axes[1])
#     axes[1].set_ylabel(r'$\Gamma_{AB}$')
#     diagonal_corrs('trace_centered', axes[2])
#     axes[2].set_ylabel(r'$\Gamma_{AB}$')
#     diagonal_corrs('time_averages', axes[3])
#     axes[3].set_ylabel(r'$\Gamma_{AB}$')
#
#     axes[0].set_title('A', x=-.15, fontsize='xx-large')
#     axes[1].set_title('B', x=-.15, fontsize='xx-large')
#     axes[2].set_title('C', x=-.15, fontsize='xx-large')
#     axes[3].set_title('D', x=-.15, fontsize='xx-large')
#     # plt.show()
#     plt.savefig('{}/Figure 3.png'.format(args.figs_location), dpi=300)
#     plt.close()


def main(args):
    """ The KL-Divergence """
    
    def kl_div_plot():
        # Here we start plotting both KL-Divergences, The good part with this is that it takes into consideration the standard deviation as well
        # sns.barplot(data=kld_df, x='combo', y='kldiv 1', hue='dataset', hue_order=['S', 'NL', 'C'], order=[r'$\phi$', r'$f$', r'$\Delta$', r'$\tau$', r'$x_0$', r'$x_\tau$', r'$\alpha$'], ci=33,
        #             ax=ax, edgecolor='black')
        sns.boxplot(data=kld_df.replace(symbols['physical_units']), x='variable', y='value', hue='dataset', hue_order=datasets, order=hierarchy, showfliers=False, palette=cmap)
        plt.legend(title='')
        plt.xlabel('')
        plt.ylabel(r'$D_{KL}$')
    
    """ diagonal pair plots """
    
    def diagonal_corrs(kind, ax):
        
        # print(over_time_and_lineages.columns)
        #
        # print(np.unique(over_time_and_lineages.AB_parameter_pair))
        #
        # # where we store what to plot
        # plot_df = pd.DataFrame(columns=over_time_and_lineages.columns)
        #
        # # # put in the diagonal parameter pair correlations
        # # for label in np.unique(over_time_and_lineages.AB_parameter_pair):
        # #     plot_df = plot_df.append(over_time_and_lineages[(over_time_and_lineages['AB_parameter_pair'] == label)], ignore_index=True)
        #
        # plot_df = plot_df[plot_df['kind'] == kind].copy()
        
        sns.barplot(data=over_time_and_lineages[over_time_and_lineages['kind'] == kind], x='AB_parameter_pair', y='correlation', hue='dataset', hue_order=datasets,
                    order=symbols[kind].values(), ci=95, dodge=True, ax=ax, edgecolor='black', palette=cmap)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.get_legend().remove()
        ax.set_ylim([-.05, .38])
        # plt.show()
        # plt.close()
    
    """ the pearson correlation coefficient between all same-generation lineage pair-bacteria in each dataset separately in Physical Units """
    
    def diagonal_units(ax, filename='physical_units', n_boots=500):
        # This is where we will save the dataframe keeping the correlations
        directory = '{}/{} correlation dataframes, {} bootstraps'.format(filename, n_boots)
        
        # the one we will plot
        corr_df = pd.DataFrame(columns=['AB_parameter_pair', 'correlation', 'dataset'])
        
        # Gather the data into a bigger dataframe
        for variable in phenotypic_variables:
            # We will add this pair parameter selection
            to_add = pd.read_csv(directory + '/{}, {}.csv'.format(variable, variable))
            print('to_add')
            print(to_add)
            to_add['AB_parameter_pair'] = r'{}'.format(symbols['physical_units'][variable], symbols['physical_units'][variable])
            
            # contains all the correlations of all pair phenotypic_variables
            corr_df = corr_df.append(to_add, ignore_index=True).reset_index(drop=True)
        
        print('corr df')
        print(corr_df)
        
        # # Stylistic purposes
        # sns.set_context('paper')
        # sns.set_style("whitegrid", {'axes.grid': True})
        
        # plot this dataframe
        # ax.yaxis.grid(which='major')
        sns.boxplot(
            x='AB_parameter_pair',
            y='correlation',
            hue='dataset',
            hue_order=datasets,
            showfliers=False,
            order=hierarchy,
            data=corr_df,
            ax=ax,
            palette=cmap
        )
        ax.axhline(0, ls='-', color='black')
        ax.set_xlabel('')
        ax.set_ylabel(r'$\rho(\cdot, \cdot)$')  # , fontsize='x-large'
        ax.set_xticklabels(ax.get_xticklabels())  # , fontsize='x-large'
        # ax.legend(title='') # , fontsize='large'
        ax.get_legend().remove()
        # plt.show()
        # plt.close()
        
    # read the csv file
    kld_df = pd.read_csv('{}/{}'.format(args.save_folder, args.pair_kld))  # 'Data/kullback_leibler_divergences_for_pair_lineages.csv'
    over_time_and_lineages = pd.read_csv('{}/{}'.format(args.save_folder, args.lin_and_time))  # 'Data/bootstrapped pearson correlations of time-averages.csv'
    
    seaborn_preamble()
    
    # create the figure and construct the layout of the figure
    cmap = sns.color_palette('tab10')[5:8]
    
    # # plot the figures
    # kl_div_plot()
    # plt.tight_layout()
    # plt.savefig('{}/pair KL divergence.png'.format(args.figs_location), dpi=300)
    # plt.close()
    #
    # fig, ax = plt.subplots(tight_layout=True)
    # diagonal_corrs('physical_units', ax)
    # ax.set_ylabel(r'$\Gamma_{AB}$')
    # plt.legend()
    # plt.savefig('{}/pair physical units gamma.png'.format(args.figs_location), dpi=300)
    # plt.close()
    #
    # fig, ax = plt.subplots(tight_layout=True)
    # diagonal_corrs('time_averages', ax)
    # ax.set_ylabel(r'$\Gamma_{AB}$')
    # plt.legend()
    # plt.savefig('{}/pair time-averages gamma.png'.format(args.figs_location), dpi=300)
    # plt.close()
    #
    # fig, ax = plt.subplots(tight_layout=True)
    # diagonal_corrs('trace_centered', ax)
    # ax.set_ylabel(r'$\Gamma_{AB}$')
    # plt.legend()
    # plt.savefig('{}/pair trace-centered gamma.png'.format(args.figs_location), dpi=300)
    # plt.close()

    fig, axes = plt.subplots(1, 3, tight_layout=True, sharey=True, figsize=[7, 4])
    axes[0].set_ylabel(r'$\Gamma_{AB}$')
    diagonal_corrs('physical_units', axes[0])
    diagonal_corrs('time_averages', axes[1])
    diagonal_corrs('trace_centered', axes[2])
    axes[0].get_legend()
    # axes[1].get_legend().remove()
    # axes[2].get_legend().remove()
    axes[0].set_title('A', x=-.15, fontsize='xx-large')
    axes[1].set_title('B', x=-.15, fontsize='xx-large')
    axes[2].set_title('C', x=-.15, fontsize='xx-large')
    plt.savefig('{}/pair gamma all three.png'.format(args.figs_location), dpi=300)
    plt.close()
    
    
if __name__ == '__main__':
    import argparse
    import os
    import time
    import numpy as np
    
    # The final_generations for the different datasets
    fg = {
        'SM': 45,
        'lambda_LB': 100,
        'Maryam_LongTraces': 70,
        'MG1655_inLB_LongTraces': 200,
        'LAC_M9': np.nan
    }
    
    # How long does running this take?
    first_time = time.time()

    data_origin = 'SM'

    print(data_origin)
    
    parser = argparse.ArgumentParser(description='Create the artificial lineages, ergodicity breaking parameters, and the KL Divergences.')
    parser.add_argument('-data_origin', '--data_origin', metavar='', type=str, help='What is the label for this data for the Data and Figures folders?', required=False, default=data_origin)
    parser.add_argument('-save', '--save_folder', metavar='', type=str, help='Where to save the dataframes.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Data/' + data_origin)
    parser.add_argument('-pair_kld', '--pair_kld', metavar='', type=str,
                        help='The filename of the dataframe that contains the kullback-leibler divergences between distributions of all pair lineages.',
                        required=False, default='kullback_leibler_divergences_for_pair_lineages.csv')
    parser.add_argument('-lat', '--lin_and_time', metavar='', type=str, help='The filename of the dataframe that contains the model between distributions of all pair lineages over lineages and time.',
                        required=False, default='over_lineages_and_time.csv')
    parser.add_argument('-f', '--figs_location', metavar='', type=str, help='Where the figures are saved.',
                        required=False, default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/Figures/' + data_origin)
    args = parser.parse_args()

    main(args)

    # How long did it take to do the mother machine?
    sm_time = time.time()

    print("--- took {} mins in total (SM data) ---".format(sm_time / 60))
