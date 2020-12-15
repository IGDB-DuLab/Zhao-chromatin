# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 19:49:38 2019

@author: ZZG
"""

import os
import copy
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter
import matplotlib.pyplot as plt
from itertools import combinations, product
from scipy.spatial.distance import pdist, squareform


plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 16
plt.rcParams['xtick.major.pad'] = 6

file_path = r'./data'

CAL_matrix = pd.read_csv(os.path.join(file_path, r'CAL.txt'), sep='\t', index_col=0)
CAL_matrix_binary = copy.deepcopy(CAL_matrix)

binary_mask = CAL_matrix > 0
CAL_matrix_binary[binary_mask] = 1
CAL_matrix_binary[~binary_mask] = 0

exp_matrix_cell = CAL_matrix_binary.sum(axis=1)
exp_matrix_strain = CAL_matrix_binary.sum(axis=0)


terminal_cell_list = list(CAL_matrix.index.values)


# Proportion of integration sites at which the GFP is expressed in a cell

fig, ax = plt.subplots(figsize=(6,4))
np.linspace(0, 100, 20)
range_list = [-0.01, 5, 10,  15, 20,  25, 30,  35, 40,  45, 50, 55,  60, 65,  70,  75, 80, 85,  90, 95,  100]
strain_number_list = [i/100*113 for i in range_list]

plt.hist(exp_matrix_cell, bins=strain_number_list, color='#52b25b', width=5.15, linewidth=1, edgecolor='lightgray')
yticks = [np.round(i/100*364, 0) for i in [3,6,9,12]]
plt.yticks(yticks, [3,6,9,12])
xticks = list([i*113/100 for i in np.arange(0, 105, 20)])
plt.xticks(xticks, [np.round(i/113, 1) for i in xticks])
plt.xlim([0,113])
plt.subplots_adjust(left=0.105, right=0.955, top=0.94 , bottom=0.1)



# Proportion of GFP expression cells at an integration site

fig, ax = plt.subplots(figsize=(6,4))
range_list = [-0.01, 5, 10,  15, 20,  25, 30,  35, 40,  45, 50, 55,  60, 65,  70,  75, 80, 85,  90, 95,  100]
cell_number_list = [i/100*364 for i in range_list]

plt.hist(exp_matrix_strain, bins=cell_number_list, color='#52b25b', width=18, linewidth=1, edgecolor='lightgray')
yticks = [np.round(i/100*113, 0) for i in [0,4,8,12,16]]
plt.yticks(yticks, [0,4,8,12,16])
xticks = list([i*360/100 for i in np.arange(0, 105, 20)])
plt.xticks(xticks, [np.round(i/360, 1) for i in xticks])
plt.xlim([-2,366])
plt.subplots_adjust(left=0.105, right=0.955, top=0.94 , bottom=0.1)




############################################################        information content        #################################################################

def information_content(category_dict, category):
    
    ic = 0
    for i in category:
        if i in category_dict.keys():
            percent_ = category_dict[i] / sum(category_dict.values())
            ic -= percent_*np.log(percent_)
        else:
            pass
    return ic

all_value = CAL_matrix.unstack()

cell_num = len(CAL_matrix.index)
position_num = len(CAL_matrix.columns)

min_value = np.nanmin(all_value)
max_value = np.nanmax(all_value)

# Between cells

simul = []
for num in [5,10,15,20,25,30,35,40,45,50]:
    
    percent_bin = num
#    percent_bin = 10
    percent_range = np.linspace(min_value,max_value, percent_bin+1)
    percent_range[0] = percent_range[0] - 0.01
#    percent_range[-1] = np.nanmax(all_value)
    category = np.arange(1,percent_bin+1)
    
#    all_value_rank = pd.qcut(all_value,  percent_bin, labels=np.arange(1,percent_bin+1))
    all_value_rank = pd.cut(all_value,  percent_range, labels=category)
    all_value_rank_matrix = all_value_rank.unstack().T
    
    max_dict_cell = {}
    for i in category:
        max_dict_cell[i] = int(cell_num/percent_bin)
    max_ic = information_content(max_dict_cell, category)
    
    information_content_all_cell = []
    
    for i in all_value_rank_matrix.columns.values:
        
        cell_ = all_value_rank_matrix[i]
        cell_value_count_ = dict(Counter(cell_.values))
        ic_ = information_content(cell_value_count_, category)
        information_content_all_cell.append([i, ic_])
        
    information_content_all_cell_pd = pd.DataFrame(information_content_all_cell)
    information_content_all_cell_pd.columns = ['cell', 'IC']
    
    information_content_all_cell_pd['IC'].max()/max_ic
    information_content_all_cell_pd['IC'].median()/max_ic
    information_content_all_cell_pd['IC'].mean()/max_ic
    
    simul.append([num, information_content_all_cell_pd['IC'].max()/max_ic, \
                  information_content_all_cell_pd['IC'].median()/max_ic, information_content_all_cell_pd['IC'].mean()/max_ic])


simul_pd = pd.DataFrame(simul)
simul_pd.columns = ['bin_num', 'max', 'median', 'mean']

fig, ax = plt.subplots(figsize=(4,4))
ax = sns.barplot(x='bin_num', y='median', data=simul_pd, color='dodgerblue')
plt.xticks(np.arange(10), np.arange(5,55, 5))
plt.ylim(0,1)
plt.ylabel('Median / Max')
plt.xlabel('Number of bins')
plt.subplots_adjust(left=0.185, right=0.975, top=0.915 , bottom=0.165)

# Between positions
simul = []
for num in [5,10,15,20,25,30,35,40,45,50]:
    
    percent_bin = num
    percent_range = np.linspace(0,max_value, percent_bin+1)
    percent_range[0] = percent_range[0] - 0.01
    percent_range[-1] = np.nanmax(all_value)   
    
    category = np.arange(1,percent_bin+1)
#    all_value_rank = pd.qcut(all_value,  percent_bin, labels=np.arange(1,percent_bin+1))
    all_value_rank = pd.cut(all_value,  percent_range, labels=category)
    all_value_rank_matrix = all_value_rank.unstack().T
    
    max_dict_cell = {}
    for i in category:
        max_dict_cell[i] = int(position_num/percent_bin)
    max_ic = information_content(max_dict_cell, category)
    
    information_content_all_cell = []
    for i in all_value_rank_matrix.index.values:
        
        cell_ = all_value_rank_matrix.loc[i]
        cell_value_count_ = dict(Counter(cell_.values))
        ic_ = information_content(cell_value_count_, category)
        information_content_all_cell.append([i, ic_])
        
    information_content_all_cell_pd = pd.DataFrame(information_content_all_cell)
    information_content_all_cell_pd.columns = ['cell', 'IC']
    
    information_content_all_cell_pd['IC'].max()/max_ic
    information_content_all_cell_pd['IC'].median()/max_ic
    information_content_all_cell_pd['IC'].mean()/max_ic
    
    simul.append([num, information_content_all_cell_pd['IC'].max()/max_ic, \
                  information_content_all_cell_pd['IC'].median()/max_ic, information_content_all_cell_pd['IC'].mean()/max_ic])


simul_pd = pd.DataFrame(simul)
simul_pd.columns = ['bin_num', 'max', 'median', 'mean']


fig, ax = plt.subplots(figsize=(4,4))
ax = sns.barplot(x='bin_num', y='median', data=simul_pd, color='dodgerblue')
plt.xticks(np.arange(10), np.arange(5,55, 5))
plt.ylim(0,1)
plt.xlabel('Number of bins')
plt.ylabel('Median / Max')
plt.subplots_adjust(left=0.185, right=0.975, top=0.915 , bottom=0.165)

    
############################################################        information content        #################################################################




exp_matrix = pd.read_csv(os.path.join(file_path, r'embryo_matrix.txt'), sep='\t', index_col=0)
exp_matrix_binary = copy.deepcopy(exp_matrix)
exp_matrix_mask = exp_matrix > 0
exp_matrix_binary[exp_matrix_mask] = 1
exp_matrix_binary[~exp_matrix_mask] = 0

divergence_matrix = pd.DataFrame(squareform(pdist(exp_matrix.T)), index = exp_matrix.columns, columns = exp_matrix.columns)


exp_matrix_binary_terminal = exp_matrix_binary.loc[terminal_cell_list]



strain_list = list(exp_matrix.columns.values)
strain_list = list(set([i.split('_')[0] for i in strain_list]))
strain_list.sort()

output_divergence = []
for each_strain in strain_list:
    
    emb_ = exp_matrix_binary_terminal.iloc[:, exp_matrix_binary_terminal.columns.str.contains(each_strain)]
    emb_exp_cell_ = emb_.mean(axis=1).sum()
    
    emb_list = list(emb_.columns.values)
    emb_combinations_ = list(combinations(emb_list, 2))
    
    for each in emb_combinations_:
        
        emb_sum_ = emb_[list(each)].sum(axis=1)
        
        nonconsistent_ = emb_sum_[emb_sum_ == 1]
        nonconsistent_ratio_ = len(nonconsistent_) / len(emb_sum_)
        
        divergence_ = divergence_matrix.loc[each[0]][each[1]]
        
        output_divergence.append([each_strain, each[0], each[1], divergence_])

output_divergence_pd = pd.DataFrame(output_divergence)
output_divergence_pd.columns = ['strain', 'emb1', 'emb2', 'euclidean']

output_divergence_pd = output_divergence_pd.set_index('strain')
output_divergence_pd = output_divergence_pd[['euclidean']]
output_corr_pd_mean = output_divergence_pd.mean(level=0)


output_other_divergence = []
for each_strain in strain_list:
    
    emb_ = exp_matrix_binary_terminal.iloc[:, exp_matrix_binary_terminal.columns.str.contains(each_strain)]
    
    other_emb = exp_matrix_binary_terminal.iloc[:, ~exp_matrix_binary_terminal.columns.str.contains(each_strain)]
    
    
    emb_list = list(emb_.columns.values)
    other_emb_list = list(other_emb.columns.values)
    emb_combinations_ = list(product(emb_list, other_emb_list))
    
    for each in emb_combinations_:
        
        emb_sum_ = exp_matrix_binary_terminal[list(each)].sum(axis=1)
        
        divergence_ = divergence_matrix.loc[each[0]][each[1]]
        
        output_other_divergence.append([each_strain, each[0], each[1], divergence_])

output_other_divergence_pd = pd.DataFrame(output_other_divergence)
output_other_divergence_pd.columns = ['strain', 'emb1', 'emb2','euclidean']

output_other_divergence_pd = output_other_divergence_pd.set_index('strain')
output_other_divergence_pd = output_other_divergence_pd[['euclidean']]
output_other_divergence_pd_mean = output_other_divergence_pd.mean(level=0)



# plot
input_other = output_other_divergence_pd_mean
input_within = output_corr_pd_mean

input_other['flag'] = 'inter'
input_within['flag'] = 'intra'

combine_sns = pd.concat([input_within, input_other])

fig, ax = plt.subplots(figsize=(4, 4))
ax=sns.histplot(combine_sns, x="euclidean", kde=True, hue="flag", edgecolor='lightgray', binwidth=2,palette=[ 'dodgerblue', 'red'] )
ax.get_legend().remove()
plt.subplots_adjust(left=0.165, right=0.95, top=0.92 , bottom=0.16)


