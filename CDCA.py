# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 21:40:23 2019

@author: ZZG
"""


import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram

# tranfroming dict to matrix type
def dict_to_matrix(dict_input):
    
    index_list = []
    columns_list = []
    
    if type(list(dict_input.values())[0]) == dict:
        
        index_list_tmp = []
        for key,values in dict_input.items():
            index_list_tmp.extend(values.keys())
            columns_list.append(key)
        index_list = list(set(index_list_tmp))
    
        matrix_list = []
        for i in columns_list:        # each row
        
            count_each_columns = []
            for j in index_list:       # each col
                if j in dict_input[i].keys():
                    count_each_columns.append(dict_input[i][j])
                else:
                    count_each_columns.append(0)
            
            matrix_list.append(count_each_columns)    
    
    elif type(list(dict_input.values())[0]) == list:
        
        index_list_tmp = []
        for key,values in dict_input.items():
            index_list_tmp.extend(values)
            columns_list.append(key)
        index_list = list(set(index_list_tmp))
        
        matrix_list = []
        for i in columns_list:        # each row
        
            count_each_columns = []
            for j in index_list:       # each col
                if j in dict_input[i]:
                    count_each_columns.append(1)
                else:
                    count_each_columns.append(0)
            
            matrix_list.append(count_each_columns)
    
    matrix_pd = pd.DataFrame(matrix_list)
    matrix_pd = matrix_pd.T
    
    matrix_pd.columns = columns_list
    matrix_pd.index = index_list

    return matrix_pd
    
# acquiring element of each cluster
def get_cluster_classes(den, label='ivl'):
    cluster_idxs = {}
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs.setdefault(c, [])
                cluster_idxs[c].append(int(i))

    cluster_classes = {}
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l

    return cluster_classes

plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 10

red_white_blue = ['#1662ad','w','red'] # red-yellow--white-blue
red_white_blue.reverse()
cmap = LinearSegmentedColormap.from_list('red_white_blue', red_white_blue[::-1], N=100)
cm.register_cmap("red_white_blue", cmap)
cpal = sns.color_palette("red_white_blue",  n_colors=100)


file_path = r'./data'

CAL_matrix = pd.read_csv(os.path.join(file_path, r'CAL.txt'), sep='\t', index_col=0)
terminal_cell_list = list(CAL_matrix.index.values)

Chromatin_co_dynamic_divergence = pd.DataFrame(squareform(pdist(CAL_matrix.T)), index=CAL_matrix.columns, columns=CAL_matrix.columns)


# identifying Chromatin co-dynamic region
color_dict = {'lime':1,'darkorange':2, 'g':3, '#0000ff':4, 'gold':5, 'red':6, 'm':7, \
              'dodgerblue':8, 'yellowgreen':9, 'grey':np.nan}
hierarchy.set_link_color_palette(['lime','darkorange', 'g','#0000ff', 'gold', 'red', 'm', 'dodgerblue', 'yellowgreen'])


link = linkage(Chromatin_co_dynamic_divergence, "average")
dendro = dendrogram(link, p=50, color_threshold=65, labels=Chromatin_co_dynamic_divergence.index, above_threshold_color='grey')

cluster_dict = get_cluster_classes(dendro)
cluster_matrix = dict_to_matrix(cluster_dict)

cluster_output = []
for color_ in cluster_matrix.columns:
    
    cluster_ = cluster_matrix[color_]
    cluster_ = cluster_[cluster_==1]
    for i in cluster_.index.values:
        cluster_output.append([i, color_dict[color_]])

cluster_output_pd = pd.DataFrame(cluster_output)
cluster_output_pd.columns = ['Chromosome location (Mb)', 'Cluster']
cluster_output_pd = cluster_output_pd.set_index('Chromosome location (Mb)') 
cluster_output_pd = cluster_output_pd.dropna()

cluster_final = cluster_output_pd.reindex(dendro['ivl'])  
cluster_final = cluster_final.replace(10, 'NA')


# plot graph

top90 = np.percentile(Chromatin_co_dynamic_divergence.unstack().dropna(), 95)
bottom10 = np.percentile(Chromatin_co_dynamic_divergence.unstack().dropna(), 5)

g = sns.clustermap(Chromatin_co_dynamic_divergence, method='average', metric='euclidean',
                  row_cluster=True, col_cluster=True,
                  linewidths=0, xticklabels=False, yticklabels=False,
                  cmap=cpal, robust=True, vmin=bottom10, vmax=top90)
