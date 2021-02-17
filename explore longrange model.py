#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names
)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.transforms as mtransforms
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os
from numpy.fft import rfft, irfft, rfftfreq
from scipy import signal, fftpack
import random
from statsmodels.tsa.seasonal import seasonal_decompose


def fourier():
    for lin_id in pu.lineage_ID.unique():
        print(lin_id)
        rel = pu[pu['lineage_ID'] == lin_id].copy().dropna().sort_values('generation')
        
        dt = 0.05
        if discrete_time:
            t = np.arange(0, np.sum(rel['generationtime']), dt)
            print(t)
        
        ran = np.random.normal(0, 1, len(rel))
        
        for variable in ['length_birth', 'ran', 'growth_rate', 'generationtime', 'division_ratio']:
            n = len(t)
            
            if variable != 'ran':
                fhat = np.fft.fft(((rel[variable] - rel[variable].mean()) / rel[variable].std()).values, n)
                # fhat = np.fft.fft(rel[variable].values, n)
            else:
                fhat = np.fft.fft((ran - np.mean(ran)) / np.std(ran), n)
            psd = fhat * np.conj(fhat) / n
            frequencies = (1 / (dt * n)) * np.arange(n)
            L = np.arange(1, np.floor(n / 2), dtype='int')
            
            print(frequencies, fhat, psd, L, sep='\n---\n')
            
            plt.plot(frequencies[L], psd[L], label=variable)
            plt.legend()
            plt.show()
            plt.close()
            
            
def longtraces():
    def croscor():
        distance_array = np.arange(0, 21)
    
        count = 0
    
        print(list(combinations(['length_birth', 'growth_rate', 'generationtime', 'division_ratio'], 2)))
    
        for var1, var2 in tuple(combinations(['length_birth', 'growth_rate', 'generationtime', 'division_ratio'], 2)):
            print(var1, var2)
            v1 = rel[var1].values
            v2 = rel[var2].values
            plt.plot(distance_array, [pearsonr(v1[:-d], v2[d:])[0] if d != 0 else pearsonr(v1, v2)[0] for d in distance_array], color=cmap[count], alpha=1, ls='-',
                     label=r'{}[:-d] {}[d:]'.format(symbols['physical_units'][var1], symbols['physical_units'][var2]))
            plt.plot(distance_array, [pearsonr(v2[:-d], v1[d:])[0] if d != 0 else pearsonr(v2, v1)[0] for d in distance_array], color=cmap[count], alpha=.5, ls='--')
            count += 1
            
    def timeseries():
        plt.axhline(0, ls='--', c='k')
        for variable in ['length_birth', 'growth_rate', 'generationtime', 'division_ratio']:

            plt.plot(rel['start_time'].values, (rel[variable] - rel[variable].mean()), label=symbols['physical_units'][variable])
        
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MG1655_inLB_LongTraces/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    # pu['fold_growth'] = np.exp(pu['fold_growth'])
    art = shuffle_info(pu, True)
    discrete_time = True
    
    for lin_id in pu.lineage_ID.unique():
        print(lin_id)
        rel = pu[pu['lineage_ID'] == lin_id].copy().dropna().sort_values('generation')
        
        print(rel.columns)

        lin = rel[['length_birth', 'growth_rate', 'generationtime', 'division_ratio']].copy()
        lin = (((lin - lin.mean()) / lin.std()) + 5).values
        result = seasonal_decompose(lin, period=10, model='multiplicative')
        
        fig, ax = plt.subplots(4, 1, figsize=[5, 5])
        
        # print(result.trend)
        # print(result.seasonal)
        # print(result.resid)
        # print(result.observed)
        result.plot()
        plt.show()
        plt.close()

        # croscor()
        # timeseries()
        # plt.legend()
        # plt.show()
        # plt.close()


longtraces()
exit()

pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/MG1655_inLB_LongTraces/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
# pu['fold_growth'] = np.exp(pu['fold_growth'])
art = shuffle_info(pu, True)
discrete_time = True

for dataset in pu.dataset.unique():
    for trap_id in pu[pu['dataset'] == dataset].trap_ID.unique():
        for trace in ['A', 'B']:
            pass
    
