#!/usr/bin/env bash

import pandas as pd
import numpy as np
from scipy.stats import linregress, zscore
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, symbols, limit_lineage_length, trace_center_a_dataframe
from string import ascii_uppercase as uppercase_letters


def main():
    # import the data
    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/physical_units.csv')

    physical_units = limit_lineage_length(physical_units, min_gens=50)
    trace_centered = trace_center_a_dataframe(physical_units)
    
    # trace_centered = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/trace_centered.csv')
    population_sampled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/population_lineages.csv')
    shuffled_pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/shuffled_generations.csv')
    tc_shuffled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/tc_shuffled.csv')
    variation_of_expanding_mean = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/variation_of_expanding_mean.csv')
    
    kinds = ['Trace', 'Population', 'Trace-Centered', 'Shuffled', 'Shuffled TC']
    dataframes = [physical_units, population_sampled, trace_centered, shuffled_pu, tc_shuffled]
    
    for variable in phenotypic_variables:
        
        fig, axes = plt.subplots(nrows=2, ncols=3, tight_layout=True, figsize=[7, 4.66])
        axes = axes.flatten()
        
        # Set the titles saying what kind of lineage we are looking at
        # and also show the mean of all the lineages across all generation-time
        for ax, df, kind in zip(axes[:-1], dataframes, kinds):
            ax.axhline(df[variable].mean(), ls='-', color='black')
            ax.set_title(kind)
            ax.set_xlim([1, 50])
            ax.yaxis.grid(True)
            ax.set_xlabel(r'$N$')
            ax.set_ylabel(symbols['time_averages'][variable]+r'$_{N}$' if variable not in ['length_birth', 'length_final'] else symbols['time_averages'][variable]+r'$_{,N}$')
    
        # plot each lineage
        for dataset in physical_units['dataset'].unique():

            for trap_id in physical_units[(physical_units['dataset'] == dataset)]['trap_ID'].unique():

                for trace in ['A', 'B']:

                    apply_mask = lambda df : df[(df['trap_ID'] == trap_id) & (df['dataset'] == dataset) & (df['trace'] == trace)][variable].copy().expanding().mean().values

                    for df, ax in zip(dataframes, axes):
                        lineage = apply_mask(df)

                        # start plotting
                        ax.plot(np.arange(1, 1 + len(lineage)), lineage, marker=',')
                  
        # show the variance analysis for the 5 different kinds of lineages
        for label, marker, color in zip(kinds, ['o', 'x', '^', '+', 'v'], sns.color_palette('tab10')):

            # the variance and generations of the TAs of this type of lineage
            dist = variation_of_expanding_mean[(variation_of_expanding_mean['param'] == variable) & (variation_of_expanding_mean['label'] == label)].sort_values('generation')

            # print(dist)
            
            # So each kind of lineage starts at the 1 and keeps their slope
            dist['var'] = dist['var'] / dist[dist['generation'] == 1]['var'].values
            
            # print(dist)
            # exit()
            
            slope, intercept, _, _, std_err = linregress(np.log(dist['generation'].values), np.log(dist['var'].values))
            axes[-1].scatter(dist['generation'], dist['var'], label=kind+r'$: {:.2} \pm {:.2}$'.format(slope, std_err), color=color, marker=marker, alpha=.7)
            axes[-1].plot(dist['generation'], np.exp([intercept + slope * np.log(gen) for gen in dist['generation']]), ls='--', color=color)

            # Add the grid
            axes[-1].grid(True)
            # ax.legend()
            axes[-1].set(xscale="log", yscale="log")
            # ax.tick_params(direction='out')
            axes[-1].set_xlabel(r'$log(N)$')
            axes[-1].set_ylabel(r'$Var($'+symbols['time_averages'][variable]+r'$_{N}) / Var($'+symbols['time_averages'][variable]+'$_{0})$' if variable not in ['length_birth', 'length_final'] else
                                r'$Var($'+symbols['time_averages'][variable]+r'$_{,N}) / Var($'+symbols['time_averages'][variable]+'$_{,0})$')
        
        plt.savefig('{}.png'.format(variable), dpi=300)
        # plt.show()
        plt.close()
        


def main1(args):
    def put_all_graphs_into_a_big_grid(df, label, variables=phenotypic_variables, remove_outliers=False, suffix='No suffix'):
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        
        print(df)
        
        for count, dataset in enumerate(df['dataset'].unique()):
            print(dataset)
            pair_a = df[(df['dataset'] == dataset) & (df['trace'] == 'A')].sort_values(['trap_ID', 'generation']).copy()[phenotypic_variables].rename(columns=latex_symbols)
            pair_b = df[(df['dataset'] == dataset) & (df['trace'] == 'B')].sort_values(['trap_ID', 'generation']).copy()[phenotypic_variables].rename(columns=latex_symbols)
            
            # Symmetrize them
            a_trace = pair_a.append(pair_b, ignore_index=True).reset_index(drop=True)
            b_trace = pair_b.append(pair_a, ignore_index=True).reset_index(drop=True)
            
            # take out the outliers
            if remove_outliers:
                before_a = len(pair_a)
                before_b = len(pair_b)
                pair_a = pair_a[(np.abs(zscore(pair_a, nan_policy='omit')) < 4).all(axis=1)].reset_index(drop=True)
                pair_b = pair_b[(np.abs(zscore(pair_b, nan_policy='omit')) < 4).all(axis=1)].reset_index(drop=True)
                after_a = len(pair_a)
                after_b = len(pair_b)
                print('percentage that stays after removing all values farther than the fourth standard deviation:', after_a / before_a, after_b / before_b)
            
            # now plot
            fig, axes = plt.subplots(nrows=len(variables) - 1, ncols=len(variables) - 1, figsize=[7, 7])
            
            for row, row_var in zip(range(axes.shape[0]), list(latex_symbols.values())[1:]):
                for col, col_var in zip(range(axes.shape[1]), list(latex_symbols.values())[:-1]):
                    
                    # define the axis we will be plotting on
                    ax = axes[row, col]
                    
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
                    if col > row:
                        ax.set_frame_on(False)
                        ax.set_xticks([])
                        ax.set_yticks([])
                    else:
                        sns.regplot(x=a_trace[col_var], y=b_trace[row_var], ax=ax, line_kws={'color': 'red'}, scatter_kws={'alpha': .1, 'color': 'grey'})
                        ax.annotate(str(np.corrcoef(a_trace[col_var], b_trace[row_var])[0, 1])[:4], xy=(.5, .72), xycoords=ax.transAxes, fontsize=13, ha='center', va='bottom', color='red',
                                    bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7))
                    
                    if row == axes.shape[0] - 1:
                        ax.set_xlabel(col_var)
                    else:
                        ax.set_xlabel('')
                        ax.set_xticklabels([])
                    if col == 0:
                        ax.set_ylabel(row_var)
                    else:
                        ax.set_ylabel('')
                        ax.set_yticklabels([])
            
            plt.tight_layout(pad=.5)
            plt.savefig('{}/{}/{} {} {}.png'.format(args.figs_location, args.pc, label, dataset, suffix), dpi=300)
            # plt.show()
            plt.close()
    
    def heatmap_analogs(df, label, variables=phenotypic_variables, annot=False, suffix='No suffix'):
        
        sns.set_context('paper')
        sns.set_style("ticks", {'axes.grid': True})
        
        latex_symbols = {variable: symbols[label][variable] for variable in variables}
        
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=[7, 2.5])
        axes = axes.flatten()
        
        for count, dataset in enumerate(df['dataset'].unique()):
            print(dataset)
            pair_a = df[(df['dataset'] == dataset) & (df['trace'] == 'A')].sort_values(['trap_ID', 'generation'])
            pair_b = df[(df['dataset'] == dataset) & (df['trace'] == 'B')].sort_values(['trap_ID', 'generation'])
            
            corr_df = pd.DataFrame(columns=variables, index=variables, dtype=float)
            
            # Symmetrize them
            a_trace = pair_a.append(pair_b, ignore_index=True).reset_index(drop=True)
            b_trace = pair_b.append(pair_a, ignore_index=True).reset_index(drop=True)
            
            repeats = []
            for var_a in variables:
                for var_b in variables:
                    if var_b not in repeats:
                        
                        corr_df[var_a].loc[var_b] = np.corrcoef(a_trace[var_a], b_trace[var_b])[0, 1]
                repeats.append(var_a)
            
            corr_df = corr_df[variables].rename(columns=latex_symbols, index=latex_symbols)
            
            mask = np.ones_like(corr_df)
            mask[np.tril_indices_from(mask)] = False
            
            if count == 2:
                cbar_ax = fig.add_axes([.91, .1, .03, .8])
                sns.heatmap(corr_df, annot=annot, square=True, vmin=-1, vmax=1, mask=mask, center=0, cbar=True, ax=axes[count], cbar_ax=cbar_ax, cbar_kws={"orientation": "vertical"}, fmt='.2f')
            else:
                sns.heatmap(corr_df, annot=annot, square=True, vmin=-1, vmax=1, mask=mask, center=0, cbar=False, ax=axes[count], fmt='.2f')
            axes[count].set_title(uppercase_letters[count], x=-.2, fontsize='xx-large')
            axes[count].set_xlabel(dataset)
        
        fig.tight_layout(rect=[0, 0, .9, 1])
        plt.savefig('{}/{}/{} {}.png'.format(args.figs_location, args.pch, label, suffix), dpi=300)
        plt.show()
        plt.close()
    
    # Create the folder where we will save the figures
    create_folder('{}/{}'.format(args.figs_location, args.pc))
    create_folder('{}/{}'.format(args.figs_location, args.pch))
    
    # import the labeled measured bacteria in trace-centered units
    physical_units = pd.read_csv('{}/{}'.format(args.save_folder, args.puc)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    trace_centered = pd.read_csv('{}/{}'.format(args.save_folder, args.tcc)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    time_averages = pd.read_csv('{}/{}'.format(args.save_folder, args.tac)).sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    unique_ta = time_averages.drop_duplicates().sort_values(['dataset', 'trap_ID', 'trace'], ascending=[1, 1, 1])
    
    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        put_all_graphs_into_a_big_grid(df, kind, variables=phenotypic_variables, remove_outliers=True, suffix='full')
    
    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        put_all_graphs_into_a_big_grid(df, kind, variables=['fold_growth', 'generationtime', 'length_birth', 'growth_rate'], remove_outliers=True, suffix='main four')
    
    exit()
    
    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        heatmap_analogs(df, kind, variables=phenotypic_variables, annot=False, suffix='full')
    
    for df, kind in zip([physical_units, trace_centered, time_averages, unique_ta], ['physical_units', 'trace_centered', 'time_averages', 'unique_ta']):
        print(kind)
        heatmap_analogs(df, kind, variables=['fold_growth', 'generationtime', 'length_birth', 'growth_rate'], annot=True, suffix='main four')
    
    exit()
    
    print('pu')
    heatmap_analogs(physical_units, 'physical_units', variables=phenotypic_variables, annot=False)
    print('tc')
    heatmap_analogs(trace_centered, 'trace_centered', variables=phenotypic_variables, annot=False)
    print('ta')
    heatmap_analogs(time_averages, 'time_averages', variables=phenotypic_variables, annot=False)
    print('ta unique')
    heatmap_analogs(unique_ta, 'unique_ta', variables=phenotypic_variables, annot=False)
    
    print('pu')
    put_all_graphs_into_a_big_grid(physical_units, 'physical_units', remove_outliers=True)
    print('tc')
    put_all_graphs_into_a_big_grid(trace_centered, 'trace_centered', remove_outliers=True)
    print('ta')
    put_all_graphs_into_a_big_grid(time_averages, 'time_averages', remove_outliers=True)
    print('unique ta')
    put_all_graphs_into_a_big_grid(unique_ta, 'unique_ta', remove_outliers=True)


if __name__ == '__main__':
    main()
