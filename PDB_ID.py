#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 19:15:24 2021

@author: qz886
"""

import os

####################### Write list from BGSU file ##################################
with open('nrlist_3.209_all.csv', 'r') as f:
    lines = f.readlines()
    
ID = []

for l in lines:
    ID.append(l.split(',')[1].split('|')[0].split('\"')[1])

ID = list(set(ID))
    
with open('list_file.txt', 'w') as f:
    for i in ID:
        f.write(i+',')
        
        
######################## Clean the PDB and CIF files ##################################        
with open('list_file.txt', 'r') as f:
    lines = f.readlines()
    
ID = lines[0].split(',')
Fail = []

for id in ID:
    PDBexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/' + id + '.pdb.gz')
    CIFexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/' + id + '.cif.gz')
    
    if PDBexist:
        if CIFexist:
            os.system('rm /Users/qz886/Desktop/Dual_Library_Update/' + id + '.cif.gz')
    
    elif CIFexist:
        continue
    
    else:
        print('No '+id)
        Fail.append(id)
        

for id in ID:
    PDBexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/' + id + '.pdb')
    CIFexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/' + id + '.cif')
    
    if PDBexist:
        os.system('rm /Users/qz886/Desktop/Dual_Library_Update/' + id + '.pdb.gz')
    
    elif CIFexist:
        os.system('rm /Users/qz886/Desktop/Dual_Library_Update/' + id + '.cif.gz')
    
    else:
        print('No '+id)
        Fail.append(id)


######################## Get equivalence class from BGSU file ##################################
EquivalenceClass = {}

with open('nrlist_3.209_all.csv', 'r') as f:
    lines = f.readlines()

for l in lines:
    repstr = l.split(',')[1].split('\"')[1].split('\"')[0]
    rep = repstr.split('|')[0] + '_'
    repstr = repstr.split('+')
    for s in repstr:
        rep += s.split('|')[2]
        
    members = l.split('\"')[5].split(',')
    mset = []
    for m in members:
        mID = m.split('|')[0] + '_'
        m = m.split('+')
        for ms in m:
            mID += ms.split('|')[2]
        mset.append(mID)
    EquivalenceClass[rep] = mset



####################### Get representitive IDs for the FSE RNA PDB files ##################################    
with open('FSE_list.txt', 'r') as f:
    lines = f.readlines()
    
listCoV = lines[0].split(', ')
repCoV = []

for ID in listCoV:
    found = False
    for rep in EquivalenceClass:
        if not found:
            for x in EquivalenceClass[rep]:
                if ID in x:
                    repCoV.append(rep.split('_')[0])
                    found = True
    if not found:
        print(ID)

# repCoV = list(set(repCoV))
        
with open('FSE_rep_list.txt', 'w') as f:
    for rep in repCoV:
        f.write(rep+'\n')



####################### Get representitive IDs for Riboswitch PDB files ##################################    
with open('Riboswitch_list.txt', 'r') as f:
    lines = f.readlines()
   
listCoV = lines[0].split(', ')
repCoV = []

for ID in listCoV:
    found = False
    for rep in EquivalenceClass:
        if not found:
            for x in EquivalenceClass[rep]:
                if ID in x and not found:
                    repCoV.append(rep.split('_')[0])
                    found = True
    if not found:
        print(ID)
# repCoV = list(set(repCoV))
       
with open('Riboswitch_rep_list.txt', 'w') as f:
    for rep in repCoV:
        f.write(rep+'\n')
        
