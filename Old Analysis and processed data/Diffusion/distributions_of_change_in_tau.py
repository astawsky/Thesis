import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import seaborn as sns
from CustomFuncsAndVars.global_variables import phenotypic_variables, symbols, create_folder


def main():
    physical_units = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/Data/physical_units.csv')

    box = np.array([])
    box1 = np.array([])

    for dataset in np.unique(physical_units['dataset']):
    
        for trap_id in np.unique(physical_units[(physical_units['dataset'] == dataset)]['trap_ID']):
        
            for trace in ['A', 'B']:
                
                lineage = physical_units[(physical_units['trap_ID'] == trap_id) & (physical_units['dataset'] == dataset) & (physical_units['trace'] == trace)].copy()['generationtime']
                
                box = np.append(box, lineage[1:].values - lineage[:-1].values)

                box1 = np.append(box1, np.abs(lineage[1:].values - lineage[:-1].values))
    
    box = np.array(box).flatten()
    box1 = np.array(box1).flatten()
    
    print(box)
    print(box1)

    sns.distplot(box)
    plt.show()
    plt.close()

    sns.distplot(box1)
    plt.show()
    plt.close()
    
    param = 'generationtime'
    
    cmap = sns.color_palette('tab10')


if __name__ == '__main__':
    main()
