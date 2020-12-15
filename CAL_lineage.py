# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:56:57 2019

@author: ZZG
"""

import os
import copy
import itertools 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations, product
from scipy.stats import t, mannwhitneyu, pearsonr
from scipy.spatial.distance import pdist, squareform
from statsmodels.sandbox.stats.multicomp import multipletests

### transforming cell lineage name
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

# plot mean + std graph
def mean_std(data_pd=None, x=None, y=None, order = None, hue=None, hue_color=None, s=10, color='dodgerblue', ylim=[0,1], legend=False, ax=None):

    list_samples=[] # making a list of arrays
    for i in order:
        list_samples.append(list(data_pd[data_pd[x] == i][y].values))
    
    def W_array(array): # function that returns W based on the array provided
        
        W = np.std(array, ddof=1)
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
        
    else:
    
        for key, values in data_pd.groupby(hue):
            
            list_samples=[] # making a list of arrays

            for i in order:
                list_samples.append(list(values[values[x] == i][y].values))   
                
            W_list = list()
            mean_list = list()
            for i in range(len(list_samples)):
                W_list.append(W_array(list_samples[i])) # makes a list of W for each array
                mean_list.append(np.nanmean(list_samples[i])) # same for the means to plot
        	
            color_ = hue_color[key]
            
            ax.errorbar(x=order, y=mean_list, yerr=W_list, fmt='o', color=color_, markersize=s, label=key)
            #plt.axvline(.25, ls='--') # this is only to demonstrate that 95%
                                      # of the 95% CI contain the actual mean

        if legend:
            plt.legend(frameon=False)
#        plt.legend(frameon=False)
        plt.xticks(order, order) 
        plt.ylim(ylim)            
        
    return ax    



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



# extracting the upper triangular matrix
def triu_matrix(matrix_input):
    
    len_matrix = len(matrix_input)
    model_triu = pd.DataFrame(np.triu(pd.DataFrame(np.eye(len_matrix)).replace(1,2).replace(0,1).replace(2,0)))
    model_triu = (model_triu == 1)
    model_triu.index = matrix_input.index
    model_triu.columns = matrix_input.columns
    matrix = matrix_input[model_triu]
    
    return matrix


# p-value correction for multiple tests
def pval_correct(pvalue_list):
    all_cell_pvalue_df = pd.DataFrame(pvalue_list)
    all_cell_pvalue_df.columns = ['pvalue']  # values

    pvalue_df = all_cell_pvalue_df[['pvalue']]
    mask = np.isfinite(pvalue_df)
    pval_corrected = np.full(pvalue_df.shape, np.nan)
    pval_corrected[mask] = multipletests(pvalue_df.values[mask], method='fdr_bh')[1]
    pval_corrected_df = pd.DataFrame(pval_corrected, index = pvalue_df.index, columns = pvalue_df.columns)
    all_cell_pvalue_df['qvalue'] = pval_corrected_df
    
    return all_cell_pvalue_df


file_path = r'./data'

# prepare data
relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")

funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}






CAL_matrix = pd.read_csv(os.path.join(file_path, r'CAL.txt'), sep='\t', index_col=0)
fate_similarity = pd.read_csv(os.path.join(file_path, r'fate_similarity_position_df.txt'), sep='\t', index_col=0)

terminal_cell_list = list(CAL_matrix.index.values)

CAL_divergence_pd = pd.DataFrame(squareform(pdist(CAL_matrix)), index=CAL_matrix.index.values, columns=CAL_matrix.index.values)


cell_combinations = list(itertools.combinations(list(CAL_matrix.index.values), 2))
cell_pair_lineage_distance_and_divergence = [] 
count = 0
for each_pair in cell_combinations:
    
    cell_1 = each_pair[0]
    cell_2 = each_pair[1]
    
    if len(cell_name_transfer(cell_1)) == len(cell_name_transfer(cell_2)):

        cell_lineage_distance_ = cell_lineage_distance(cell_1, cell_2)
        cell_exp_distance = CAL_divergence_pd.loc[cell_1, cell_2]
        
        fate_divergence_ = 1 - fate_similarity.loc[cell_1, cell_2]
        cell_pair_lineage_distance_and_divergence.append([cell_1, cell_2, cell_lineage_distance_, cell_exp_distance, fate_divergence_])
        
    count += 1
    print(count)
        
cell_pair_lineage_distance_and_divergence_pd = pd.DataFrame(cell_pair_lineage_distance_and_divergence)
cell_pair_lineage_distance_and_divergence_pd.columns = ['cell_1', 'cell_2', 'Cell lineage distance', 'Chromatin landscape divergence', 'Fate divergence']    

cell_pair_lineage_distance_and_divergence_pd = cell_pair_lineage_distance_and_divergence_pd[cell_pair_lineage_distance_and_divergence_pd['Cell lineage distance'] % 2 == 0]

#cell_pair_lineage_distance_and_divergence_pd['Chromatin landscape divergence'].mean()


# plot graph
plt.rcParams['xtick.major.pad'] = 2
plt.rcParams["errorbar.capsize"] = 5
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 16

fig, ax = plt.subplots(figsize=(5,3.3))

fate_pd = copy.deepcopy(cell_pair_lineage_distance_and_divergence_pd)
fate_pd['Cell lineage distance'] = fate_pd['Cell lineage distance'] * 3
fate_lineage_distance_order = list(set(fate_pd['Cell lineage distance'].values))
fate_lineage_distance_order.sort()
mean_std(fate_pd, x='Cell lineage distance', y="Chromatin landscape divergence", order=fate_lineage_distance_order,\
         s=9, color='dodgerblue', ylim=[-2.8,27], ax = ax)
plt.yticks([0,5,10,15,20,25], ['0','5','10','15','20','25'])

ax2=ax.twinx()

lineage_pd = copy.deepcopy(cell_pair_lineage_distance_and_divergence_pd)
lineage_pd['Cell lineage distance'] = lineage_pd['Cell lineage distance']*3 + 2
lineage_distance_order = list(set(lineage_pd['Cell lineage distance'].values))
lineage_distance_order.sort()
ax = mean_std(lineage_pd, x='Cell lineage distance', y="Fate divergence", order=lineage_distance_order,\
         s=9, color='red', ylim=[-0.09,1.39], ax = ax2)
plt.subplots_adjust(left=0.120, right=0.865, top=0.895 , bottom=0.15)
plt.yticks([0.0,0.25,0.5,0.75,1.0,1.25], ['0.00','0.25','0.50','0.75','1.00','1.25'])
x = [i-1 for i in lineage_distance_order]
plt.xticks(x, ['2','4','6','8','10','12','14','16','18'])




# CAL bin
CAL_bin = [0, 10, 20, 1000]

cell_pair_lineage_distance_and_divergence_pd['bin'] = pd.cut(cell_pair_lineage_distance_and_divergence_pd['Chromatin landscape divergence'], bins=CAL_bin, labels=[1,2,3])

plt.rcParams["errorbar.capsize"] = 10
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 16

lineage_distance_bin = list(set(cell_pair_lineage_distance_and_divergence_pd['Cell lineage distance']))
lineage_distance_bin.sort()

for each_cld in [6,8,10,12,14]:
    
    fate_and_lineage_pd_ = cell_pair_lineage_distance_and_divergence_pd[cell_pair_lineage_distance_and_divergence_pd['Cell lineage distance'] == each_cld]
    
    CAL_bin_list_ = list(set(fate_and_lineage_pd_['bin']))
    CAL_bin_list_.sort()
    
    final_bin = []
    for i in CAL_bin_list_:
        
        number_ = len(fate_and_lineage_pd_[fate_and_lineage_pd_['bin'] == i])
        if number_ > 3:
            final_bin.append(i)
        
    mean_IC(x='bin', y='Fate divergence', data_pd=fate_and_lineage_pd_, color='dodgerblue', \
                 order=CAL_bin_list_, IC=0.95, ylim=[0.2,1], s=15, figsize=(4,4))
    
    plt.subplots_adjust(left=0.155, right=0.92, top=0.880 , bottom=0.110)
    plt.xlim(0.5,3.5)    


# CAL landscape transient points
CAL_divergence_pd_normalize_cld = pd.read_csv(os.path.join(file_path, r'CAL_normalize_cld.txt'), sep='\t', index_col=0)

terminal_cell_list_replace = cell_name_transfer(list(CAL_divergence_pd_normalize_cld.index.values))

CAL_divergence_pd_normalize_cld.index = terminal_cell_list_replace
CAL_divergence_pd_normalize_cld.columns = terminal_cell_list_replace

cell_list = pd.read_csv(os.path.join(file_path, r'cell-list.txt'), sep='\t', index_col=0)
cell_list = list(cell_list.index.values)
all_cell_list_replace = cell_name_transfer(cell_list)
all_cell_list_replace = copy.deepcopy(all_cell_list_replace)
all_cell_list_replace.append('L')
all_cell_list_replace.append('L1')
all_cell_list_replace.append('L0')

fate_similarity_change_index = copy.deepcopy(fate_similarity)
fate_similarity_change_index.index = cell_name_transfer(list(fate_similarity_change_index.index))
fate_similarity_change_index.columns = cell_name_transfer(list(fate_similarity_change_index.columns))


cutoff = 0.05
fc_cutoff = 2
number_cutoff = 3

each_site_strain = {}

all_cell_stats_pvalue_count = []
for cell in all_cell_list_replace:
    
    child_1 = cell+'0'
    child_2 = cell+'1'
    
    child_1_list = []
    child_2_list = []
    mother_list = []
    
    for i in terminal_cell_list_replace:
        
        if child_1 in i:
            child_1_list.append(i)
        if child_2 in i:
            child_2_list.append(i)
        if cell in i:
            mother_list.append(i)
            
    if len(child_1_list) < number_cutoff or len(child_2_list) < number_cutoff:
        
        count = np.nan
        compare_count = np.nan
        sig_strain = np.nan
        strain_compare = np.nan
        diff_fate = np.nan
        
    else:
        
        diff_fate = 1-fate_similarity_change_index.loc[child_1][child_2]
        
        intra1_ = list(combinations(child_1_list, 2))
        intra2_ = list(combinations(child_2_list, 2))
        inter_ = list(product(child_1_list, child_2_list))
        
        intra1_distance_ = []
        for i in intra1_:
            intra1_distance_.append(CAL_divergence_pd_normalize_cld.loc[i[0]][i[1]])
            
        intra2_distance_ = []
        for i in intra2_:
            intra2_distance_.append(CAL_divergence_pd_normalize_cld.loc[i[0]][i[1]])    
            
        inter_distance_ = []
        for i in inter_:
            inter_distance_.append(CAL_divergence_pd_normalize_cld.loc[i[0]][i[1]])
    
        intra1_mean = np.nanmean(intra1_distance_)
        intra2_mean = np.nanmean(intra2_distance_)            
        inter_mean = np.nanmean(inter_distance_)            

        H1, pvalue1 = mannwhitneyu(intra1_distance_, inter_distance_, alternative='two-sided')
        H2, pvalue2 = mannwhitneyu(intra2_distance_, inter_distance_, alternative='two-sided')

        all_cell_stats_pvalue_count.append([cell_name_transfer(cell), diff_fate, intra1_mean, intra2_mean, inter_mean, pvalue1, pvalue2])              # values
    
all_cell_stats_pvalue_df = pd.DataFrame(all_cell_stats_pvalue_count)

all_cell_stats_pvalue_df.columns = ['cell_name', 'Fate divergence','intra1_mean','intra2_mean', 'inter_mean', 'pvalue1', 'pvalue2']  # values
all_cell_stats_pvalue_df['Transition score of CAL'] = (all_cell_stats_pvalue_df['inter_mean']/all_cell_stats_pvalue_df['intra1_mean'] + \
                                            all_cell_stats_pvalue_df['inter_mean']/all_cell_stats_pvalue_df['intra2_mean'])/2


all_cell_stats_pvalue_df['qvalue1'] = pval_correct(all_cell_stats_pvalue_df['pvalue1'])['qvalue']
all_cell_stats_pvalue_df['qvalue2'] = pval_correct(all_cell_stats_pvalue_df['pvalue2'])['qvalue']

all_cell_stats_pvalue_df['left'] = np.where(((all_cell_stats_pvalue_df['intra1_mean'] < all_cell_stats_pvalue_df['inter_mean']) & (all_cell_stats_pvalue_df['qvalue1'] < cutoff)), 1,0)
all_cell_stats_pvalue_df['right'] = np.where(((all_cell_stats_pvalue_df['intra2_mean'] < all_cell_stats_pvalue_df['inter_mean']) & (all_cell_stats_pvalue_df['qvalue2'] < cutoff)), 1,0)

all_cell_stats_pvalue_df['flag'] = all_cell_stats_pvalue_df['left'] + all_cell_stats_pvalue_df['right']

sig_cell = list(all_cell_stats_pvalue_df[all_cell_stats_pvalue_df['flag']!=0]['cell_name'])
other_cell = list(set(all_cell_stats_pvalue_df['cell_name']) -set(sig_cell))

len(sig_cell)

# plot graph
pearsonr(all_cell_stats_pvalue_df['Fate divergence'], all_cell_stats_pvalue_df['Transition score of CAL'])
sns.regplot(x='Transition score of CAL', y='Fate divergence',data = all_cell_stats_pvalue_df, color='dodgerblue')
plt.xlim(0.65,1.65)
plt.xticks([0.8,1.0,1.2,1.4,1.6], ['0.8','1.0','1.2','1.4','1.6'])
plt.yticks([0.0,0.2,0.4,0.6,0.8,1.0], ['0.0','0.2','0.4','0.6','0.8','1.0'])
plt.subplots_adjust(left=0.225, right=0.900, top=0.965 , bottom=0.150)



