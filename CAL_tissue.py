# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 11:02:37 2019

@author: ZZG
"""

import os
import copy
import numpy as np
import pandas as pd
from scipy.stats import t
from collections import Counter
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from scipy.spatial.distance import pdist, squareform

# transforming cell lineage name
def except_cell(cell_name):
    
    transform_name = None
    
    if cell_name[0] == 'L':
        
        len_cell_name = len(cell_name)
        
        for i in range(1,len_cell_name+1):
            
            if cell_name[:-i] in list(relpace_name_to_cell_name_sheet.index.values):
                
                transform_name = (relpace_name_to_cell_name_sheet.loc[[cell_name[:-i]]].values[0])[0]
                break
            
        append_name = cell_name[len_cell_name-i:]
        
        for j in append_name:
            
            if j == '0':
                transform_name+='a'
            else:
                transform_name+='p'
                
    else:
        
        len_cell_name = len(cell_name)
        
        for i in range(1,len_cell_name+1):
            
            if cell_name[:-i] in list(cell_name_to_relpace_name_sheet.index.values):
                
                transform_name = (cell_name_to_relpace_name_sheet.loc[[cell_name[:-i]]].values[0])[0]
                break
            
        append_name = cell_name[len_cell_name-i:]
        
        for j in append_name:
            
            transform_name += transform[j]
    
    return transform_name
        



def cell_name_transfer(cell_list):
    
    transfer_cell_list = []
    if len(cell_list) == 0:
        pass
        
    elif type(cell_list) != list:
        if cell_list[0] == "L":
            
            if cell_list in list(relpace_name_to_cell_name_sheet.index.values):
            
                transfer_cell_list = (relpace_name_to_cell_name_sheet.loc[[cell_list]].values[0])[0]
                
            else:
                
                transfer_cell_list = except_cell(cell_list)
        else:
            if cell_list in list(cell_name_to_relpace_name_sheet.index.values):
                
                transfer_cell_list = (cell_name_to_relpace_name_sheet.loc[[cell_list]].values[0])[0]
            
            else:
                transfer_cell_list = except_cell(cell_list)

                
            
    else:
        if cell_list[0][0] == "L":
            for i in cell_list:
                
                if i in list(relpace_name_to_cell_name_sheet.index.values):
            
                    transfer_cell_list.append((relpace_name_to_cell_name_sheet.loc[[i]].values[0])[0])
                
                else:
                
                    transfer_cell_list.append(except_cell(i))

            
        else:
            for i in cell_list:
                
                if i in list(cell_name_to_relpace_name_sheet.index.values):
                    
                    transfer_cell_list.append((cell_name_to_relpace_name_sheet.loc[[i]].values[0])[0])
                    
                else:
                    
                    transfer_cell_list.append(except_cell(i))

    return transfer_cell_list


# calculating cell lineage distance
def cell_lineage_distance(cell_1, cell_2):
    
    if cell_1[0] == 'L':
        pass
    else:
        cell_1 = cell_name_transfer(cell_1)
        cell_2 = cell_name_transfer(cell_2)
    
    
    
    cell_1_list = list(cell_1)
    cell_2_list = list(cell_2)
    
    len_cell_1_list = len(cell_1_list)
    len_cell_2_list = len(cell_2_list)
    
    min_len = min(len_cell_1_list, len_cell_2_list)
    
    lowest_common_ancestor=''
    for i in range(min_len):
        if cell_1_list[i] == cell_2_list[i]:
            lowest_common_ancestor += cell_1_list[i]
        else:
            break
            
    distance_1 = len_cell_1_list - len(lowest_common_ancestor)
    distance_2 = len_cell_2_list - len(lowest_common_ancestor)
        
    distance = distance_1+distance_2
    
    return distance


# plot mean + CI graph
def mean_IC(data_pd=None, x=None, y=None, order = None, hue=None, hue_color=None, IC=0.95, s=10, color='dodgerblue', ylim=[0,1], figsize=(8,6), legend=False):

    list_samples=[] # making a list of arrays
#    lineage_distance = list(set(data_pd[x]))
#    lineage_distance.sort()
    for i in order:
        list_samples.append(list(data_pd[data_pd[x] == i][y].values))
    
    def W_array(array, conf=IC): # function that returns W based on the array provided
        t_num = t(df = len(array) - 1).ppf((1 + conf) /2)
        W = t_num * np.std(array, ddof=1) / np.sqrt(len(array))
        return W # the error
    
    if hue == None:
    
        W_list = list()
        mean_list = list()
        for i in range(len(list_samples)):
            W_list.append(W_array(list_samples[i])) # makes a list of W for each array
            mean_list.append(np.nanmean(list_samples[i])) # same for the means to plot
    		
        fig, ax = plt.subplots(figsize=figsize)
        
        ax.errorbar(x=order, y=mean_list, yerr=W_list, fmt='o', color=color, markersize=s)
        #plt.axvline(.25, ls='--') # this is only to demonstrate that 95%
                                  # of the 95% CI contain the actual mean
        plt.xticks(order, order) 
        plt.ylim(ylim)
        
    else:
        
        fig, ax = plt.subplots(figsize=figsize)
        for key, values in data_pd.groupby(hue):
            
            list_samples=[] # making a list of arrays

            for i in order:
                list_samples.append(list(values[values[x] == i][y].values))   
                
            W_list = list()
            mean_list = list()
            for i in range(len(list_samples)):
                W_list.append(W_array(list_samples[i], IC)) # makes a list of W for each array
                mean_list.append(np.nanmean(list_samples[i])) # same for the means to plot
        	
            color_ = hue_color[key]
            
            ax.errorbar(x=order, y=mean_list, yerr=W_list, fmt='o', color=color_, markersize=s, label=key)
            #plt.axvline(.25, ls='--') # this is only to demonstrate that 95%
                                      # of the 95% CI contain the actual mean
        if legend:
            plt.legend(frameon=False)
        plt.xticks(order, order) 
        plt.ylim(ylim)            
        
    return ax  


def mean_IC_add(ax = None, data_pd=None, x=None, y=None, order = None, hue=None, hue_color=None, IC=0.95, s=10, color='dodgerblue', ylim=[0,1], figsize=(8,6), legend=False):

    
    list_samples=[] # making a list of arrays
#    lineage_distance = list(set(data_pd[x]))
#    lineage_distance.sort()
    for i in order:
        list_samples.append(list(data_pd[data_pd[x] == i][y].values))
    
    def W_array(array, conf=IC): # function that returns W based on the array provided
        t_num = t(df = len(array) - 1).ppf((1 + conf) /2)
        W = t_num * np.std(array, ddof=1) / np.sqrt(len(array))
        return W # the error
    
    if hue == None:
    
        W_list = list()
        mean_list = list()
        for i in range(len(list_samples)):
            W_list.append(W_array(list_samples[i])) # makes a list of W for each array
            mean_list.append(np.nanmean(list_samples[i])) # same for the means to plot

        ax.errorbar(x=order, y=mean_list, yerr=W_list, fmt='o', color=color, markersize=s)
        #plt.axvline(.25, ls='--') # this is only to demonstrate that 95%
                                  # of the 95% CI contain the actual mean
        plt.xticks(order, order) 
        plt.ylim(ylim)
        
    return ax  


# extracting the upper triangular matrix
def triu(matrix):
    
    matrix_bk = copy.deepcopy(matrix)
    mask = np.zeros_like(matrix, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    mask_pd = pd.DataFrame(mask)
    mask_pd.index = matrix.index
    mask_pd.columns = matrix.columns
    matrix_bk = matrix[mask_pd]
    
    return matrix_bk

plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'arial'
plt.rcParams["errorbar.capsize"] = 10


file_path = r'./data'



# prepare data

relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")

funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}



CAL_matrix = pd.read_csv(os.path.join(file_path, r'CAL.txt'), sep='\t', index_col=0)
terminal_cell_list = list(CAL_matrix.index.values)

CAL_divergence_pd = pd.DataFrame(squareform(pdist(CAL_matrix)), index=CAL_matrix.index.values, columns=CAL_matrix.index.values)



cld_matrix = pd.read_csv(os.path.join(file_path, r'ALL_CLD_MATRIX.txt'), sep='\t', index_col=0)
cld_matrix = cld_matrix.loc[terminal_cell_list][terminal_cell_list]
cld_matrix = cld_matrix.astype(float)


np.fill_diagonal(CAL_divergence_pd.values, np.nan)
np.fill_diagonal(cld_matrix.values, np.nan)



#fate dict
fate = pd.read_csv(os.path.join(file_path, r'cell-fate.txt'), sep='\t', index_col=0)
fate = fate['fate']
fate = fate.dropna()


fate_list = ['Neu', 'Pha', 'Ski', 'Mus', 'Int']
fate_cell_list = {}
for each_fate in fate_list:
    tmp_ = list(fate[fate == each_fate].index.values)
    fate_cell_list[each_fate] = tmp_




all_pair_same_fate_cls_pd = []
for each_fate in fate_list:
    
    fate_cell_list_ = fate_cell_list[each_fate]
    distance_ = CAL_divergence_pd.loc[fate_cell_list_][fate_cell_list_]
    cld_ = cld_matrix.loc[fate_cell_list_][fate_cell_list_]
    
    distance_ = triu(distance_)
    cld_ = triu(cld_)
    
    combine = pd.concat([distance_.unstack().dropna(), cld_.unstack().dropna()], axis=1)
    combine['fate_flag'] = each_fate
    
    if type(all_pair_same_fate_cls_pd) == list:
        all_pair_same_fate_cls_pd = combine
    else:
        all_pair_same_fate_cls_pd = pd.concat([all_pair_same_fate_cls_pd, combine])
    

all_pair_same_fate_cls_pd = all_pair_same_fate_cls_pd.dropna()
all_pair_same_fate_cls_pd.columns = ['Chromatin landscape divergence', 'Cell lineage distance', 'Tissue']


same_cld_list = list(set(list(all_pair_same_fate_cls_pd['Cell lineage distance'].values)))
all_pair_same_fate_cls_pd.reset_index(inplace=True)
all_pair_same_fate_cls_pd = all_pair_same_fate_cls_pd.set_index(['Cell lineage distance', 'Tissue'])
all_pair_same_fate_cls_pd_mean = all_pair_same_fate_cls_pd.mean(level=('Cell lineage distance', 'Tissue'))  # mean

all_pair_same_fate_cls_pd_mean.reset_index(inplace=True)
all_pair_same_fate_cls_pd_mean['Legend'] = 'Same tissue'


all_pair_same_fate_cls_pd.reset_index(inplace=True)



all_pair_same_fate_cls_pd_copy = copy.deepcopy(all_pair_same_fate_cls_pd)
all_pair_same_fate_cls_pd_copy = all_pair_same_fate_cls_pd_copy[['level_0', 'level_1']]
all_pair_same_fate_cls_pd_copy = all_pair_same_fate_cls_pd_copy.reset_index(drop=True)
all_pair_same_fate_cls_pd_copy['flag'] =1
all_pair_same_fate_cls_pd_copy = all_pair_same_fate_cls_pd_copy.set_index(['level_0', 'level_1'])
all_pair_same_fate_matrix = all_pair_same_fate_cls_pd_copy['flag'].unstack()
all_pair_same_fate_matrix = all_pair_same_fate_matrix.reindex(terminal_cell_list)
all_pair_same_fate_matrix = all_pair_same_fate_matrix.T.reindex(terminal_cell_list)
all_pair_same_fate_matrix  = all_pair_same_fate_matrix.fillna(0)

all_pair_same_fate_matrix = (all_pair_same_fate_matrix + all_pair_same_fate_matrix.T) /2
all_pair_same_fate_matrix = all_pair_same_fate_matrix.replace(0.5,1)


same_fate_flag = all_pair_same_fate_matrix == 1
diff_fate_flag = all_pair_same_fate_matrix == 0


# get same/distinct tissue cell-cell CAL divergence

cell_distance_pd_diff = CAL_divergence_pd[diff_fate_flag]
cld_matrix_diff = cld_matrix[diff_fate_flag]


cell_distance_pd_same = CAL_divergence_pd[same_fate_flag]
cld_matrix_same = cld_matrix[same_fate_flag]








########################################################             raw           ######################################


all_fate_cell_list = list(fate[~((fate == 'other') | (fate == 'Dea'))].index.values)

output_tissue = []
for fate_, cell_list_ in fate_cell_list.items():
    
    
    other_fate = list(set(all_fate_cell_list) - set(cell_list_))
    
    for each_cell in cell_list_:

        cell_distance_diff_list_ = list(cell_distance_pd_diff.loc[each_cell][other_fate].dropna().values)
        cell_distance_same_list_ = list(cell_distance_pd_same.loc[each_cell][cell_list_].dropna().values)
        
        
        output_tissue.append([fate_, each_cell, np.nanmean(cell_distance_same_list_), np.nanmean(cell_distance_diff_list_)])

output_tissue_pd = pd.DataFrame(output_tissue)
output_tissue_pd.columns = ['Tissue', 'cell', 'same', 'diff']



# plot
same_f_tissue = output_tissue_pd[['Tissue', 'cell', 'same']]
same_f_tissue['flag'] = 'same tissue'
same_f_tissue.columns = ['Tissue', 'cell', 'Chromatin landscape divergence', 'Legend']
diff_f_tissue = output_tissue_pd[['Tissue', 'cell', 'diff']]
diff_f_tissue['flag'] = 'diff tissue'
diff_f_tissue.columns = ['Tissue', 'cell', 'Chromatin landscape divergence', 'Legend']

sns_pd_tissue = pd.concat([same_f_tissue, diff_f_tissue])


sns_pd_tissue = sns_pd_tissue.replace('Neu', 1)
sns_pd_tissue = sns_pd_tissue.replace('Pha', 2)
sns_pd_tissue = sns_pd_tissue.replace('Ski', 3)
sns_pd_tissue = sns_pd_tissue.replace('Mus', 4)
sns_pd_tissue = sns_pd_tissue.replace('Int', 5)
sns_pd_tissue = sns_pd_tissue[sns_pd_tissue['Tissue'].isin([1,2,3,4,5])]


sns_pd_same = sns_pd_tissue[sns_pd_tissue['Legend'] == 'same tissue']
sns_pd_diff = sns_pd_tissue[sns_pd_tissue['Legend'] == 'diff tissue']


sns_pd_same = sns_pd_same[['Chromatin landscape divergence', 'Tissue']]
sns_pd_diff = sns_pd_diff[['Chromatin landscape divergence', 'Tissue']]

sns_pd_diff['Tissue'] = sns_pd_diff['Tissue'].astype(int)
sns_pd_diff['Tissue'] = sns_pd_diff['Tissue'] + 0.4
sns_pd_same['flag'] = 'same tissue'
sns_pd_diff['flag'] = 'diff tissue'

combine_two = pd.concat([sns_pd_same, sns_pd_diff])
new_order = list(set(combine_two['Tissue']))
new_order.sort()



color_dict = {1:'red', 2:'gold', 3:'limegreen', 4:'dodgerblue', 5:'fuchsia'}
ax = mean_IC(x='Tissue', y='Chromatin landscape divergence', hue='Tissue', hue_color= color_dict, data_pd=sns_pd_same,  order=new_order, IC=0.95, \
             ylim=[7,30], s=15, figsize=(5,5))

mean_IC_add(ax=ax,x='Tissue', y='Chromatin landscape divergence', data_pd=sns_pd_diff, color='darkgrey',  order=new_order, IC=0.95, \
             ylim=[7,30], s=15, figsize=(5,5))
plt.subplots_adjust(left=0.155, right=0.92, top=0.880 , bottom=0.110)

plt.xlim([0.7,5.7])
plt.xticks([1.2, 2.2, 3.2, 4.2, 5.2],['Neu', 'Pha', 'Ski', 'Mus', 'Int'])


for each_tissue in fate_list:

    u_, p_ = mannwhitneyu(output_tissue_pd[output_tissue_pd['Tissue'] == each_tissue]['same'], \
                          output_tissue_pd[output_tissue_pd['Tissue'] == each_tissue]['diff'], alternative='two-sided')
    print(p_)


########################################################             raw           ######################################











########################################################       remove symmetry     ######################################

symmetry_matrix = pd.read_csv(os.path.join(file_path, r'symmetry-matrix.txt'), sep='\t', index_col=0)

all_fate_cell_list = list(fate[~((fate == 'other') | (fate == 'Dea'))].index.values)

output_remove_symmetry = []
for fate_, cell_list_ in fate_cell_list.items():
    
    
    other_fate = list(set(all_fate_cell_list) - set(cell_list_))
    for each_cell in cell_list_:
        
        symmetric_cell_list = []
        symmetry__ = symmetry_matrix.loc[each_cell].dropna()
        
        if len(symmetry__) >= 1:
            symmetric_cell_list = list(symmetry__.index.values)
            
        cell_list_copy = copy.deepcopy(cell_list_)
        cell_list_copy = set(cell_list_copy) - set(symmetric_cell_list)
        
        
        cell_cld_same_list_ = list(set(list(cld_matrix_same.loc[each_cell][cell_list_copy].dropna().values)))
        combine_cld = cell_cld_same_list_


        symmetry_mask = symmetry_matrix.loc[each_cell]
        symmetry_mask = symmetry_mask.fillna(0)
        symmetry_mask = symmetry_mask.astype(bool)
        del_symmetry_mask = ~symmetry_mask
        cell_cld_same_mask_ = (cld_matrix_same.loc[each_cell].isin(combine_cld)) & (del_symmetry_mask)
        
        
        cell_distance_same_pd_ = cell_distance_pd_same.loc[each_cell][cell_cld_same_mask_].dropna()
        cld_same_pd_ = cld_matrix_same.loc[each_cell][cell_cld_same_mask_].dropna()
        same_pd_ = pd.concat([cell_distance_same_pd_, cld_same_pd_], axis=1)
        same_pd_.columns = ['distance', 'cld']
        same_pd_ = same_pd_.set_index('cld')
        same_pd_mean_ = same_pd_.mean(level=0)
        same_pd_mean_mean_ = same_pd_mean_.mean().values[0]
        
        symmetry_mask = copy.deepcopy(symmetry_matrix)
        symmetry_mask = symmetry_mask.fillna(0)
        symmetry_mask = symmetry_mask.astype(bool)
        del_symmetry_mask = ~symmetry_mask
        cell_cld_diff_mask_ = (cld_matrix_diff.isin(combine_cld)) & (del_symmetry_mask)

        
        cell_distance_diff_pd_ = cell_distance_pd_diff[cell_cld_diff_mask_].unstack().dropna()
        cld_diff_pd_ = cld_matrix_diff[cell_cld_diff_mask_].unstack().dropna()
        diff_pd_ = pd.concat([cell_distance_diff_pd_, cld_diff_pd_], axis=1)
        diff_pd_.columns = ['distance', 'cld']
        diff_pd_ = diff_pd_.set_index('cld')
        diff_pd_mean_ = diff_pd_.mean(level=0)
        diff_pd_mean_mean_ = diff_pd_mean_.mean().values[0]
        
        
        output_remove_symmetry.append([fate_, each_cell, same_pd_mean_mean_, diff_pd_mean_mean_])

output_remove_symmetry_pd = pd.DataFrame(output_remove_symmetry)
output_remove_symmetry_pd.columns = ['Tissue', 'cell', 'same', 'diff']

same_f_remove_symmetry = output_remove_symmetry_pd[['Tissue', 'cell', 'same']]
same_f_remove_symmetry['flag'] = 'same tissue'
same_f_remove_symmetry.columns = ['Tissue', 'cell', 'Chromatin landscape divergence', 'Legend']
diff_f_remove_symmetry = output_remove_symmetry_pd[['Tissue', 'cell', 'diff']]
diff_f_remove_symmetry['flag'] = 'diff tissue'
diff_f_remove_symmetry.columns = ['Tissue', 'cell', 'Chromatin landscape divergence', 'Legend']

sns_pd_remove_symmetry = pd.concat([same_f_remove_symmetry, diff_f_remove_symmetry])





# plot

sns_pd_remove_symmetry = sns_pd_remove_symmetry.replace('Neu', 1)
sns_pd_remove_symmetry = sns_pd_remove_symmetry.replace('Pha', 2)
sns_pd_remove_symmetry = sns_pd_remove_symmetry.replace('Ski', 3)
sns_pd_remove_symmetry = sns_pd_remove_symmetry.replace('Mus', 4)
sns_pd_remove_symmetry = sns_pd_remove_symmetry.replace('Int', 5)
sns_pd_remove_symmetry = sns_pd_remove_symmetry[sns_pd_remove_symmetry['Tissue'].isin([1,2,3,4,5])]


sns_pd_same_remove_symmetry = sns_pd_remove_symmetry[sns_pd_remove_symmetry['Legend'] == 'same tissue']
sns_pd_diff_remove_symmetry = sns_pd_remove_symmetry[sns_pd_remove_symmetry['Legend'] == 'diff tissue']


sns_pd_same_remove_symmetry = sns_pd_same_remove_symmetry[['Chromatin landscape divergence', 'Tissue']]
sns_pd_diff_remove_symmetry = sns_pd_diff_remove_symmetry[['Chromatin landscape divergence', 'Tissue']]

sns_pd_diff_remove_symmetry['Tissue'] = sns_pd_diff_remove_symmetry['Tissue'] + 0.4
sns_pd_same_remove_symmetry['flag'] = 'same tissue'
sns_pd_diff_remove_symmetry['flag'] = 'diff tissue'

combine_two_remove_symmetry = pd.concat([sns_pd_same_remove_symmetry, sns_pd_diff_remove_symmetry])




color_dict = {1:'red', 2:'gold', 3:'limegreen', 4:'dodgerblue', 5:'fuchsia'}
ylim = [7,18]
ax = mean_IC(x='Tissue', y='Chromatin landscape divergence', hue='Tissue', hue_color= color_dict, data_pd=sns_pd_same_remove_symmetry,  order=new_order, IC=0.95, \
             ylim=ylim, s=15, figsize=(5,5))


mean_IC_add(ax=ax,x='Tissue', y='Chromatin landscape divergence', data_pd=sns_pd_diff_remove_symmetry, color='darkgrey',  order=new_order, IC=0.95, \
             ylim=ylim, s=15, figsize=(5,5))
plt.subplots_adjust(left=0.155, right=0.92, top=0.880 , bottom=0.110)
plt.xlim([0.7,5.7])



for each_tissue in fate_list:

    u_, p_ = mannwhitneyu(output_remove_symmetry_pd[output_remove_symmetry_pd['Tissue'] == each_tissue]['same'], \
                          output_remove_symmetry_pd[output_remove_symmetry_pd['Tissue'] == each_tissue]['diff'], alternative='two-sided')
    print(p_)

########################################################       remove symmetry     ######################################









########################################################       heterogeneity     ######################################


fate_color_dict = {'Neu':'red', 'Pha':'gold', 'Ski':'limegreen', 'Mus':'dodgerblue', 'Int':'fuchsia'}


for fate_name in fate_list:

    
    fate_cell_list_ = fate_cell_list[fate_name]
    CAL_divergence_pd_ = CAL_divergence_pd.loc[fate_cell_list_][fate_cell_list_]

    
    cell_combinations = list(combinations(fate_cell_list_, 2))
    cell_pair_lineage_distance_and_exp_distance = [] 
    for each_pair in cell_combinations:
        
        cell_1 = each_pair[0]
        cell_2 = each_pair[1]
        
        if len(cell_name_transfer(cell_1)) == len(cell_name_transfer(cell_2)):
        
            
            if cell_1[0] == cell_2[0]:
                flag = cell_1[0]
            else:
                flag = np.nan

            
            cell_lineage_distance_ = cell_lineage_distance(cell_1, cell_2)
            cell_exp_distance = CAL_divergence_pd_.loc[cell_1, cell_2]
            cell_pair_lineage_distance_and_exp_distance.append([cell_1, cell_2, cell_lineage_distance_, cell_exp_distance])
            
    cell_pair_lineage_distance_and_exp_distance_df = pd.DataFrame(cell_pair_lineage_distance_and_exp_distance)
    cell_pair_lineage_distance_and_exp_distance_df.columns = ['cell_1', 'cell_2', 'Cell lineage distance', 'Chromatin landscape divergence']    
    
    cell_pair_lineage_distance_and_exp_distance_df = cell_pair_lineage_distance_and_exp_distance_df[cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'] % 2 == 0]
    
    distance_pool = list(set(list(cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'].values)))
    
    distance_pool_filter = []
    for each in distance_pool:
        
        if len(cell_pair_lineage_distance_and_exp_distance_df[cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'] == each]) >= 3:
            
            distance_pool_filter.append(each)
            
    cell_pair_lineage_distance_and_exp_distance_df = cell_pair_lineage_distance_and_exp_distance_df[cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'].isin(distance_pool_filter)]
    
    
    

    
    cld_dict = dict(Counter(cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance']))
    
    order_list = []
    for cld_, count_ in cld_dict.items():
        
        if count_ >= 3:
            order_list.append(cld_)
#    
    order_list.sort()
    
#    color_dict = [rgb2hex(int(i[0]*255), int(i[1]*255), int(i[2]*255)) for i in cpal]
    
    color_fate_ = fate_color_dict[fate_name]

    tmp = copy.deepcopy(cell_pair_lineage_distance_and_exp_distance_df)
    tmp = tmp.set_index('Cell lineage distance')
    tmp_mean = tmp['Chromatin landscape divergence'].mean(level=0)
    max_index = tmp_mean.idxmax()
    min_index = tmp_mean.idxmin()
    
    two_ = cell_pair_lineage_distance_and_exp_distance_df[cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'] == min_index]['Chromatin landscape divergence']
    max_ = cell_pair_lineage_distance_and_exp_distance_df[cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'] == max_index]['Chromatin landscape divergence']
    
    min_mean_ = two_.mean()
    min_std_ = two_.std(ddof=1)
    max_mean_ = max_.mean()
    max_std_ = max_.std(ddof=1)
    
    ymin = min_mean_ - min_std_ - 2
    ymax = max_mean_ + max_std_ + 2
    

    mean_IC(data_pd=cell_pair_lineage_distance_and_exp_distance_df, x='Cell lineage distance', y='Chromatin landscape divergence', order = order_list,\
                 hue=None, hue_color=None, IC=0.95, s=15, color=color_fate_, ylim=[0,20], figsize=(5,5), legend=False)
#    plt.xlim([1, 17])

########################################################       heterogeneity     ######################################

#
#
#order_list = [2,4,6,8,10,12,14,16]
#cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'] = cell_pair_lineage_distance_and_exp_distance_df['Cell lineage distance'].replace(18, 16)