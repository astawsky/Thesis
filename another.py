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


# # Asymptotic Approximations
# x = np.linspace(0.0001, 0.02)
# y = np.linspace(50, 400)
# plt.plot(x, 0.73*x**-.94)
# plt.plot(x, 14*x**-.45)
# plt.plot(x, np.log(2) / x, ls='--')
# plt.ylim([50, 400])
# plt.show()
# plt.close()
#
# exit()

total = pd.DataFrame()
ax = plt.axes(projection='3d')

for count, ds in enumerate([sm_datasets[-1]] + cgsc_6300_wang_exps + tanouchi_datasets + mm_datasets):
    print(ds)
    
    # pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/physical_units_without_outliers.csv')
    pu = pd.read_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
    # pu = get_time_averages_df(pu, phenotypic_variables).drop_duplicates()
    # pu.to_csv(f'/Users/alestawsky/PycharmProjects/Thesis/Datasets/{ds}/ProcessedData/z_score_under_3/time_average_without_outliers.csv')
    pu['experiment'] = ds
    pu['division_ratio'] = pu['division_ratio']
    total = total.append(pu, ignore_index=True)
    
    pu = pu.sample(n=min(100, len(pu)), replace=False)
    
    # Data for a three-dimensional line
    zline = pu['growth_rate'].values
    xline = pu['generationtime'].values
    yline = pu['division_ratio'].values
    ax.scatter3D(xline, yline, zline, alpha=.5)  # , c=cmap[count]

domain = np.linspace(- 3 * total['growth_rate'].std() + total['growth_rate'].mean(), + 3 * total['growth_rate'].std() + total['growth_rate'].mean())
fake_domain = np.array([np.log(2) / z for z in domain])

X, Z = np.meshgrid(
    fake_domain,
    domain
)

# print(X)
# print('()')
# print(Z)
# exit()

Y = 1 / np.exp(X * Z)
# Y = - (X*Z)

# surf = ax.plot_surface(X, Y, Z, alpha=.8)
y = np.array([1/np.exp(gr * gt) for gr, gt in zip(domain, fake_domain)])
# y = np.array([-(gr * gt) for gr, gt in zip(domain, fake_domain)])
ax.plot3D(fake_domain, y, domain, color='black')
ax.set_xlabel(r'$\tau$')
# ax.set_ylabel(r'f')
ax.set_ylabel(r'$f$')
# plt.yscale('log')
# ax.yaxis.set_scale('log')
ax.set_zlabel(r'$\alpha$')
plt.tight_layout()
plt.show()
plt.close()

total=total.sample(n=1000,replace=False)

print(len(total))

pcorr = np.round(pearsonr(total.dropna().generationtime.values, total.dropna().growth_rate.values)[0], 2)
print(pcorr)
print(total.columns)
sns.kdeplot(
    data=total, x='generationtime', y='growth_rate',
    label=f'{pcorr}'
)
s, i = np.round(linregress(total.generationtime.values, total.growth_rate.values)[:2], 2)
plt.plot(
    total.generationtime.unique(), [i+s*g for g in total.generationtime.unique()],
    label=f'y={i}+{s}*x'
)
plt.show()
plt.close()

# sns.regplot(
#     data=total, x='generationtime', y='growth_rate',
#     label=f'{np.round(pearsonr(total.generationtime.values, total.growth_rate.values)[0], 2)}, {np.round(linregress(total.generationtime.values, total.growth_rate.values)[0], 2)}'
# )
