# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 11:33:23 2019

@author: ZZG
"""

import os
import copy
import pandas as pd
import numpy as np


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


### sorting cell name list according lineage order
def lineage_order(cell_list):
    
    if cell_list[0][0] == "L":
        cell_list = cell_name_transfer(cell_list)
    else:
        pass
    
    cell_list_replace = cell_name_transfer(cell_list)
    cell_list_replace.sort()
    cell_list_replace = sorted(cell_list_replace, reverse=False, key=len)
    cell_list = cell_name_transfer(cell_list_replace)
    
    AB = []
    E = []
    MS = []
    D = []
    P = []
    C = []
    
    for i in cell_list:
        if i[0] == "A":
            AB.append(i)
        elif i[0] == "C":
            C.append(i)
        elif i[0] == "D":
            D.append(i)
        elif i[0] == "E":
            E.append(i)
        elif i[0] == "M":
            MS.append(i)  
        else:
            P.append(i)
    AB.sort()
    E.sort()
    MS.sort()
    D.sort()
    P.sort()
    C.sort()    
    lineage_cell_list = AB+MS+E+C+D+P  
    
    return lineage_cell_list
    

###  selecting terminall cell list
def get_terminal_cell(cell_list):
    
    if cell_list[0][0] == "L":
        cell_list = cell_name_transfer(cell_list)
    else:
        pass
    
    unique_cell_set = set(cell_list)
    unique_cell_list = list(unique_cell_set)
    transfer_cell_list = cell_name_transfer(unique_cell_list)
    transfer_cell_set = set(transfer_cell_list)
    
    not_terminal_list = []
    for i in transfer_cell_list:
        len_i = len(i)
        for j in transfer_cell_list:
            if len(j) <= len_i:
                pass
            else:
                if j[:len_i] == i:
                    not_terminal_list.append(i)
                    break
    
    
    terminal_cell_set = transfer_cell_set - set(not_terminal_list)

    terminal_cell_list = list(terminal_cell_set)
    terminal_cell_list = cell_name_transfer(terminal_cell_list)
    terminal_cell_list = lineage_order(terminal_cell_list)
    return terminal_cell_list


### compensating fluorescent intensity
def compensation(emb):
    
    embryo_name = os.path.basename(emb)
    path_ = os.path.dirname(emb)
    
    embryo_ = pd.read_csv(emb, sep='\t')
    
    cell_list = list(set(list(embryo_['cell_name'].values)))
    terminal_cell_list = get_terminal_cell(cell_list)
    middle_cell_list = list(set(cell_list) - set(terminal_cell_list))
    
    
    
    # del first and last time point
    new_emb_terminal_ = []
    for each in terminal_cell_list:
        tmp_ = embryo_[embryo_['cell_name'] == each]
        tmp_ = tmp_.sort_values(by='time')
        tmp_ = tmp_.iloc[1:,:]
        if type(new_emb_terminal_) == list:
            new_emb_terminal_ = tmp_
        else:
            new_emb_terminal_ = pd.concat([new_emb_terminal_, tmp_])
    
    new_emb_middle_ = []
    for each in middle_cell_list:
        tmp_ = embryo_[embryo_['cell_name'] == each]
        tmp_ = tmp_.sort_values(by='time')
        tmp_ = tmp_.iloc[1:-1,:]
        if type(new_emb_middle_) == list:
            new_emb_middle_ = tmp_
        else:
            new_emb_middle_ = pd.concat([new_emb_middle_, tmp_])
            
    new_emb = pd.concat([new_emb_middle_, new_emb_terminal_])
    
    
    
    # eliminate mothe cell accumulation
    g = new_emb.groupby('cell_name')
                    
    mother_ = g.nth([-1, -2])['raw-expression']
    mother_ = mother_.mean(level=0)
    
    new_emb = new_emb.set_index('cell_name')
    
    new_ = []
    for each_cell in cell_list:
        
        mother_cell_ = cell_name_transfer(cell_name_transfer(each_cell)[:-1])
        
        if mother_cell_ in cell_list:
            
            try:
                terminal_time = new_emb.loc[each_cell]['raw-expression']
#                       terminal_time = terminal_.loc[each_cell]['mean_gray_value']

                if type(terminal_time) != pd.core.series.Series:
                    new_.append([each_cell, np.nan])
                else:
                    terminal_time = terminal_time - mother_.loc[mother_cell_]
                    new_.append([each_cell, terminal_time.mean(skipna=True)])
            except:
            
                new_.append([each_cell, np.nan])

    new_pd_ = pd.DataFrame(new_)
    new_pd_.columns = ['cell', embryo_name]
    new_pd_ = new_pd_.set_index(['cell'])
    
    cutoff_strain_bin_cell = {}

    cutoff_strain_bin_cell.setdefault(embryo_name, {})
    strain_embryos_mask = copy.deepcopy(new_pd_)
    
    
    bin0_mask = (new_pd_ < 7000)
    bin1_mask = (new_pd_ >= 7000)
    
    try:
        strain_embryos_mask[bin0_mask] = 'bin0'                
    except(TypeError):
        cutoff_strain_bin_cell[embryo_name].setdefault('bin0', [])
        cutoff_strain_bin_cell[embryo_name]['bin0'] = []


    try:
        strain_embryos_mask[bin1_mask] = 'bin1'
    except(TypeError):
        cutoff_strain_bin_cell[embryo_name].setdefault('bin1', [])
        cutoff_strain_bin_cell[embryo_name]['bin1'] = []  
        


    strain_embryos_mask = strain_embryos_mask.fillna('bin0')
    
    strain_embryos_mask_max = strain_embryos_mask.max(axis=1)

    bin0_cell_tmp = list(strain_embryos_mask_max[strain_embryos_mask_max == 'bin0'].dropna().index.values)
    bin1_cell_tmp = list(strain_embryos_mask_max[strain_embryos_mask_max == 'bin1'].dropna().index.values)
    
    cutoff_strain_bin_cell[embryo_name].setdefault('bin0', [])
    cutoff_strain_bin_cell[embryo_name].setdefault('bin1', [])
    
    cutoff_strain_bin_cell[embryo_name]['bin0'].extend(bin0_cell_tmp)
    cutoff_strain_bin_cell[embryo_name]['bin1'].extend(bin1_cell_tmp)

    
    compensation_parameter = {'bin0':1, 'bin1':1.054}
    


    new_ = []
    for each_cell in cell_list:
        
        mother_cell_ = cell_name_transfer(cell_name_transfer(each_cell)[:-1])
        
        if mother_cell_ in cell_list:
            
            try:
                terminal_time = new_emb.loc[[each_cell], :]
#                       terminal_time = terminal_.loc[each_cell]['mean_gray_value']
                tmp_ = terminal_time[['raw-expression']] - mother_.loc[mother_cell_]
                tmp_.columns = ['expression']
                
                terminal_time = pd.concat([terminal_time, tmp_], axis=1)


            except:
            
                pass
            
            if type(new_) == list:
                if type(terminal_time) == pd.core.series.Series:
                    pass
                else:
                    new_ = terminal_time
            else:
                if type(terminal_time) == pd.core.series.Series:
                    pass
                else:
                    new_ = pd.concat([new_, terminal_time])
            
        
    new_.reset_index(inplace=True)
    new_ = new_.set_index('cell_name')
        

    compensation_each_embryo = []
    for each_bin in compensation_parameter.keys():
        
        bin_cell_ = cutoff_strain_bin_cell[embryo_name][each_bin]
        
        embryo_bin_ = new_.loc[new_.index.isin(bin_cell_)]

        # delta Z    
        delta_z = np.array([15.0]*len(embryo_bin_)) - embryo_bin_["Z"].values
        exp_ = embryo_bin_[["expression"]]
        
        #######################################
        negative_mask = exp_<=0
        positive_mask = exp_>0
        
        exp_positive_ =  exp_[positive_mask]["expression"] / np.power(np.array([compensation_parameter[each_bin]]*len(exp_)), delta_z)
        exp_negative_ = exp_[negative_mask]['expression']
        
        exp_positive_ = exp_positive_.fillna(0)
        exp_negative_ = exp_negative_.fillna(0)
        
        exp_compensation_ = exp_positive_ + exp_negative_
        
        exp_compensation_ = pd.DataFrame(exp_compensation_)
        exp_compensation_.columns = ['compensation-value']
        
        embryo_bin_ = pd.concat([embryo_bin_, exp_compensation_], axis=1)
        
        
        if type(compensation_each_embryo) == list:
            compensation_each_embryo = embryo_bin_
        else:
            compensation_each_embryo = pd.concat([compensation_each_embryo, embryo_bin_])


    compensation_each_embryo.reset_index(inplace=True)
    compensation_each_embryo = compensation_each_embryo.dropna()
    compensation_each_embryo = compensation_each_embryo[['cell_name', 'time', 'Z', 'raw-expression','compensation-value']]
    compensation_each_embryo.to_csv(os.path.join(path_, embryo_name+'_compensation.txt'), sep='\t', index=None)
    

    
    
    print('compensation result output succeeded!')
    
    
file_path = r'./data'

# prepare data
relpace_name_sheet = pd.read_csv(os.path.join(file_path, 'binary_sheet.txt'), sep="\t")

cell_name_to_relpace_name_sheet = relpace_name_sheet.set_index("cell_name")
relpace_name_to_cell_name_sheet = relpace_name_sheet.set_index("replace_name")

funder_cell = {'P0':'L', 'AB':'L0', 'P1':'L1', 'EMS':'L10', 'P2':'L11', 'MS':'L100', 'E':'L101', 'C':'L110', 'D':'L1110', 'P4':'L1111', 'Z2':'L11110', 'Z3':'L11111'}
transform = {'a':'0', 'p':'1', 'l':'0', 'r':'1', 'd':'0', 'v':'1'}



compensation(os.path.join(file_path, r'test_embryo.txt'))
