#!/usr/bin/env bash

from AnalysisCode.global_variables import (
    symbols, units, create_folder, phenotypic_variables, shuffle_lineage_generations, cmap, seaborn_preamble, sm_datasets,
    wang_datasets, get_time_averages_df, trace_center_a_dataframe, cgsc_6300_wang_exps, shuffle_info, lexA3_wang_exps, mm_datasets, tanouchi_datasets, dataset_names, shuffle_lineage_generations
)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress, pearsonr
from itertools import combinations
import os
from matplotlib import cm


def scatter_plot_3d(df, ax):
    # Data for a three-dimensional line
    ax.scatter3D(df['generationtime'].values, df['division_ratio'].values, df['growth_rate'].values, alpha=1)
    
    
def plot_manifold(df, ax):
    # gr_domain = np.linspace(- 3 * df['growth_rate'].std() + df['growth_rate'].mean(), + 3 * df['growth_rate'].std() + df['growth_rate'].mean())  # -3.5 *
    gt_domain = np.linspace(df['generationtime'].min(), df['generationtime'].max())
    gr_domain = np.linspace(df['growth_rate'].min(), df['growth_rate'].max())
    # gt_domain = np.linspace(- 3 * df['generationtime'].std() + df['generationtime'].mean(), + 3 * df['generationtime'].std() + df['generationtime'].mean())
    
    X, Z = np.meshgrid(
        gt_domain,
        gr_domain
    )

    Y = 1 / np.exp(X * Z)  # In the case of ln(f)
    
    ax.plot_wireframe(X, Y, Z, cmap=cm.hot, alpha=.5)  # plot_surface , alpha=.5
    
    
def plot_perfect_line(df, ax):
    # gr_domain = np.linspace(- 3 * df['growth_rate'].std() + df['growth_rate'].mean(), + 3 * df['growth_rate'].std() + df['growth_rate'].mean())
    gr_domain = np.linspace(df['growth_rate'].min(), df['growth_rate'].max())

    gt_domain = np.array([np.log(2) / z for z in gr_domain])

    y = np.array([1 / np.exp(gr * gt) for gr, gt in zip(gr_domain, gt_domain)])
    ax.plot3D(gt_domain, y, gr_domain, color='black', zorder=10000)
    
    
def script(points_are='pu'):  # can be ta
    
    total = pd.DataFrame()
    ax = plt.axes(projection='3d')
    # plt.hold(True)
    
    for ds in dataset_names[:-1]:
        if ds == 'Pooled_SM':
            continue
            
        print(ds)
        
        pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
        # pu = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
        # pu.to_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
        pu['experiment'] = ds
        # pu['division_ratio'] = np.log(pu['division_ratio'])
        
        if points_are == 'pu':
            total = total.append(pu, ignore_index=True)
            # pu = pu.sample(n=min(100, len(pu)), replace=False)
            # scatter_plot_3d(pu, ax)
        else:
            # ta = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
            try:
                ta = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
            except:
                continue
            total = total.append(ta, ignore_index=True)
            # ta = ta.sample(n=min(100, len(pu)), replace=False)
            # scatter_plot_3d(ta, ax)
        
    plot_manifold(total, ax)
    plot_perfect_line(total, ax)
    
    for ds in dataset_names[:-1]:
        if ds == 'Pooled_SM':
            continue

        pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
        # pu = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
        # pu.to_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
        pu['experiment'] = ds
        # pu['division_ratio'] = np.log(pu['division_ratio'])
        
        if points_are == 'pu':
            # total = total.append(pu, ignore_index=True)
            pu = pu.sample(n=min(500, len(pu)), replace=False)
            scatter_plot_3d(pu, ax)
        else:
            # ta = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
            try:
                ta = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
            except:
                continue
            # total = total.append(ta, ignore_index=True)
            ta = ta.sample(n=min(100, len(pu)), replace=False)
            scatter_plot_3d(ta, ax)
    
    ax.set_xlabel(r'$\tau$')
    # ax.set_ylabel(r'$\ln(f)$')
    ax.set_ylabel(r'$f$')
    ax.set_zlabel(r'$\alpha$')
    plt.tight_layout()
    plt.show()
    plt.close()
    
    
script(points_are='ta')
