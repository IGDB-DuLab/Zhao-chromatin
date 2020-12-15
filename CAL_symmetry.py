# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:28:53 2020

@author: ZZG
"""
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import t, wilcoxon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm


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


blue_white_red = ['#0070c0', (1, 1, 1), (1,0,0)] # blue-white-red2
cmap = LinearSegmentedColormap.from_list('blue_white_red', blue_white_red, N=100)
cmap.set_bad('k')
cm.register_cmap("blue_white_red", cmap)
cpal = sns.color_palette("blue_white_red",  n_colors=100)


add_color = (0,0,0)

new_cpal = [add_color]
new_cpal.extend(cpal)


plt.rcParams['font.size'] = 20
plt.rcParams['font.family'] = 'arial'
plt.rcParams["errorbar.capsize"] = 10

file_path = r'I:\position-effect\FINAL\MSB\analysis\software\scCAL\data'

# prepare data
symmetry = pd.read_csv(os.path.join(file_path, r'symmetry.txt'), sep='\t')


# symmetry predetermination

ax = sns.clustermap(symmetry[['CAL divergence', 'Mean CAL divergence between all intra-tissue cells at the same lineage distance']].T,  \
                    col_cluster=False, row_cluster=False,  xticklabels=False, yticklabels=False, \
               row_colors=None, col_colors=None, linewidth=0.01, cmap=cmap, vmin=8, vmax=15, figsize=(10, 2), linecolor='gray')

symmetry_stats = symmetry[['CAL divergence', 'Mean CAL divergence between all intra-tissue cells at the same lineage distance']]
symmetry_stats = symmetry_stats.dropna()
wilcoxon(symmetry_stats['CAL divergence'], symmetry_stats['Mean CAL divergence between all intra-tissue cells at the same lineage distance'])




# heterogeneity

cld_list = list(set(symmetry['Cell lineage distance']))

cld_list_filter = []
for each in cld_list:
    if len(symmetry[symmetry['Cell lineage distance'] == each]) >= 3:
        
        cld_list_filter.append(each)
symmetry = symmetry[symmetry['Cell lineage distance'].isin(cld_list_filter)]

cld_list_filter.sort()

mean_IC(data_pd=symmetry, x='Cell lineage distance', y='CAL divergence', order = cld_list_filter,\
             hue=None, hue_color=None, IC=0.95, s=15, color='dodgerblue', ylim=[0,20], figsize=(5,5), legend=False)
#    plt.xlim([1, 17])




