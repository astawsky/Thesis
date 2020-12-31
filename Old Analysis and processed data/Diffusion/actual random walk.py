import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
# import statsmodels.api as sm
import seaborn as sns
from scipy.stats import linregress
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder


def main():
    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/physical_units.csv')

    for dataset in np.unique(physical_units['dataset']):
    
        for trap_id in np.unique(physical_units[(physical_units['dataset'] == dataset)]['trap_ID']):
        
            for trace in ['A', 'B']:
                pu_lineage = physical_units[(physical_units['trap_ID'] == trap_id) & (physical_units['dataset'] == dataset) & (physical_units['trace'] == trace)].copy().values
                shuffled = pu_lineage.copy().sample(frac=1, replace=False).reset_index(drop=True)
                
                pu_random_walk = pu_lineage.values - pu_lineage.values[0]
                shuffled_random_walk = shuffled.values - shuffled.values[0]
                
                

    # trace_centered = pd.read_csv('Data/trace_centered.csv')
    v_cum_sum = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/variation_of_cumulative_sum.csv')
    # cum_sum_df = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/cumulative_sum.csv')
    
    where_to_save_figs = '/Users/alestawsky/PycharmProjects/Thesis/Diffusion/Figures'
    create_folder(where_to_save_figs)
    
    cmap = sns.color_palette('tab10')
    
    markers = ['.', '*', '^', '+', 'x']
    
    for param in phenotypic_variables:
        
        for count, kind in enumerate(['Population', 'Trace-Centered', 'Shuffled TC']):  # 'Shuffled', 'Trace',
            to_plot = v_cum_sum[(v_cum_sum['label'] == kind) & (v_cum_sum['param'] == param)].sort_values('generation')
            
            slope, intercept, _, _, std_err = linregress(np.log(to_plot['generation'].values), np.log(to_plot['var'].values))
            
            plt.scatter(to_plot['generation'], to_plot['var'], label=kind + r'$: {:.2} \pm {:.2}$'.format(slope, std_err), color=cmap[count], marker=markers[count], alpha=.7)
            plt.plot(to_plot['generation'], np.exp([intercept + slope * np.log(gen) for gen in to_plot['generation']]), ls='--', color=cmap[count])
        
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.title(symbols['physical_units'][param])
        plt.ylabel(r'$log(\sigma^2)$')
        plt.xlabel('log generation')
        plt.tight_layout()
        plt.savefig(where_to_save_figs + '/{}.png'.format(param), dpi=300)
        # plt.show()
        plt.close()


def the_walks():
    # import/create the trace-centered lineages
    trace_centered = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/trace_centered.csv')
    shuffled_tc = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/tc_shuffled.csv')
    population_sampled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/population_lineages.csv')
    shuffled = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/shuffled_generations.csv')
    
    param = 'generationtime'
    
    for dataset in np.unique(trace_centered['dataset']):
        
        for trap_id in np.unique(trace_centered[(trace_centered['dataset'] == dataset)]['trap_ID']):
            
            for trace in ['A', 'B']:
                tc_lineage = trace_centered[(trace_centered['trap_ID'] == trap_id) & (trace_centered['dataset'] == dataset) & (trace_centered['trace'] == trace)].copy()[param].values
                tc_shuffled_lineage = shuffled_tc[(shuffled_tc['trap_ID'] == trap_id) & (shuffled_tc['dataset'] == dataset) & (shuffled_tc['trace'] == trace)].copy()[param].values
                pop_lineage = population_sampled[(population_sampled['trap_ID'] == trap_id) & (population_sampled['dataset'] == dataset) & (population_sampled['trace'] == trace)].copy()[param].values
                shuffled_lineage = shuffled[(shuffled['trap_ID'] == trap_id) & (shuffled['dataset'] == dataset) & (shuffled['trace'] == trace)].copy()[param].values
                
                print(tc_lineage, tc_shuffled_lineage, pop_lineage, sep='\n' * 2)
                
                # plt.plot(tc_lineage, label='TC', marker='.')
                # plt.plot(tc_shuffled_lineage, label='TC Shuffled', marker='.')
                plt.plot((pop_lineage - pop_lineage.mean()).cumsum(), label='Pop', marker='.')
                plt.plot((shuffled_lineage - shuffled_lineage.mean()).cumsum(), label='shuffled', marker='.')
                plt.axhline(0, ls='--', color='black')
                
                plt.tight_layout()
                plt.legend()
                plt.show()
                plt.close()


if __name__ == '__main__':
    main()
    # the_walks()
