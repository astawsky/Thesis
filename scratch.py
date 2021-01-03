import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from CustomFuncsAndVars.global_variables import phenotypic_variables, create_folder, dataset_names



# pu = pd.read_csv('/Users/alestawsky/PycharmProjects/Thesis/ProcessedData/lambda_LB/physical_units.csv')
# param1 = 'generationtime'
# param2 = 'growth_rate'
# ta1 = [np.mean(pu[pu['lineage_ID'] == lin_id][param1]) for lin_id in pu.lineage_ID.unique()]
# ta2 = [np.mean(pu[pu['lineage_ID'] == lin_id][param2]) for lin_id in pu.lineage_ID.unique()]
# plt.scatter(ta1, ta2)
# plt.xlabel(param1)
# plt.ylabel(param2)
# plt.show()
# plt.close()

# date_mapping = {
#     '072818': '072818_SL_NL',
#     '101218': '101218_SL_NL',
#     '071318': '071318_SL',
#     '062718': '062718_SL',
#     '1012': '101218_SL_NL',
#     '1015': '1015_NL',
#     '0728': '072818_SL_NL'
# }

# for date in np.unique(list(date_mapping.values())):
#     create_folder('/Users/alestawsky/PycharmProjects/Thesis/RawData/' + date)
    
    
new_date_mapping = date_mapping.copy()


# for pair_type in ['NL', 'SL']:
#
#     files = glob.glob(r'/Users/alestawsky/PycharmProjects/Thesis/RawData/RawData/*.xls'.format(pair_type))
#     print('files', files)
#
#     for file in files:
#         print(file)
#         date = file.split('/')[-1].split('.')[0].split('_')[0]
#         new_file = r'/Users/alestawsky/PycharmProjects/Thesis/RawData/'+date_mapping[date]+'/'+file.split('/')[-1]
#
#         os.rename(file, new_file)

# sister_files = []
# neighbor_files = []
#
# for folder in np.unique(list(date_mapping.values()))[1:]:
#     print('*'*200)
#     print(folder)
#
#     files = glob.glob(r'/Users/alestawsky/PycharmProjects/Thesis/RawData/' + folder + '/*')
#
#     print(files)
#
#     for file in files:
#         if ('POS' in file) or ('Pos' in file):
#             print(file)
#             # print(file, file.split('.')[0]+'_SL'+file.split('.')[-1])
#             # os.rename(file, file.split('.')[0]+'_SL'+file.split('.')[-1])
#         else:
#             # print(file, file.split('.')[0]+'_NL'+file.split('.')[-1])
#             pass
#
#     # exit()
#     # print([file if ('POS' in file) or ('Pos' in file) else '' for file in files])
#     #
#     # print(len(files), len([file if ('POS' in file) or ('Pos' in file) else '' for file in files]))
#     # exit()

