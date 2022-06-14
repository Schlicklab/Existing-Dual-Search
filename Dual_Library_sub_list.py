#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 19:59:28 2021

@author: qz886
"""

from igraph import *
from ClassesFunctions import *
from dualGraphs import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os,sys
import re
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors


def graph_Collector(id, chains, DualGraph_lib):
    fname = id+'_'
    chainList = []
    for x in chains:
        fname += x
        chainList.append(x.split('-')[0])
    chainList = list(set(chainList))
                
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        graphs = []
        for l in lines:
            if '_' in l:
                if l != 'No subgraph for 1_1\n':
                    graphs.append(l.split()[0])
                
        graphs = list(set(graphs))
        for g in graphs:
            if g in DualGraph_lib:
                if isDNA(id, chains):
                    DualGraph_lib[g].append(fname+'*')
                else:
                    if len(chainList)>1:
                        DualGraph_lib[g].append(fname+'#')
                    else:
                        DualGraph_lib[g].append(fname)
            else:
                if isDNA(id, chains):
                    DualGraph_lib[g] = [fname+'*']
                else:
                    if len(chainList)>1:
                        DualGraph_lib[g] = [fname+'#']
                    else:
                        DualGraph_lib[g] = [fname]

    else:
        print('Cannot identify '+fname+'.txt')


# check if the chains contain DNA
def isDNA(id, chains):
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i][0:6] == 'Chains':
            if lines[i].split()[1] in chains:
                if 'DNA' in lines[i]:
                    return True
            i += 1
            while lines[i][0] == '\t':
                if lines[i].split()[0] in chains:
                    if 'DNA' in lines[i]:
                        return True
                i += 1
        else:
            i += 1
    return False


def graph_Collector_title(id, chains, DualGraph_lib):
    fname = id+'_'
    for x in chains:
        fname += x
    
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
    title = lines[0].split('\t')[3].split('\n')[0]
            
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        graphs = []
        for l in lines:
            if '_' in l:
                if l != 'No subgraph for 1_1\n':
                    graphs.append(l.split()[0])
                
        graphs = list(set(graphs))
        for g in graphs:
            if g in DualGraph_lib:
                DualGraph_lib[g].append(id)
            else:
                DualGraph_lib[g] = [id]

    else:
        print('Cannot identify '+fname+'.txt')
        
        

def graph_Collector_largePDB(id, chains, DualGraph_lib):
    fname = id+'_'
    chainList = []
    for x in chains:
        fname += x
        chainList.append(x.split('-')[0])
    chainList = list(set(chainList))
                
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        graphs = []
        for l in lines:
            if '_' in l:
                if l != 'No subgraph for 1_1\n':
                    graphs.append(l.split()[0])
                
        graphs = list(set(graphs))
        for g in graphs:
            if g in DualGraph_lib:
                if isDNA_largePDB(id, chains):
                    DualGraph_lib[g].append(fname+'*')
                else:
                    if len(chainList)>1:
                        DualGraph_lib[g].append(fname+'#')
                    else:
                        DualGraph_lib[g].append(fname)
            else:
                if isDNA_largePDB(id, chains):
                    DualGraph_lib[g] = [fname+'*']
                else:
                    if len(chainList)>1:
                        DualGraph_lib[g] = [fname+'#']
                    else:
                        DualGraph_lib[g] = [fname]

    else:
        print('Cannot identify '+fname+'.txt')


# check if the chains contain DNA
def isDNA_largePDB(id, chains):
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i][0:6] == 'Chains':
            if lines[i].split()[1] in chains:
                if 'DNA' in lines[i]:
                    return True
            i += 1
            while lines[i][0] == '\t':
                if lines[i].split()[0] in chains:
                    if 'DNA' in lines[i]:
                        return True
                i += 1
        else:
            i += 1
    return False


def graph_Collector_largePDB_title(id, chains, DualGraph_lib):
    fname = id+'_'
    for x in chains:
        fname += x
    
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
    title = lines[0].split('\t')[3].split('\n')[0]
            
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        graphs = []
        for l in lines:
            if '_' in l:
                if l != 'No subgraph for 1_1\n':
                    graphs.append(l.split()[0])
                
        graphs = list(set(graphs))
        for g in graphs:
            if g in DualGraph_lib:
                DualGraph_lib[g].append(id)
            else:
                DualGraph_lib[g] = [id]

    else:
        print('Cannot identify '+fname+'.txt')
        
        

########################### Get graph library for all structures #################################
with open('/Users/qz886/Desktop/Dual_Library_Update/list_file.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')

# These two pdb files have unusual structures, no vertex
manual = ['7BPG', '7BPF']
for ma in manual:
    ID.remove(ma)

largePDB = []
largePDBfrag = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i][0:17] == 'Too many helices:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            largePDBfrag.append(s)
            
    if lines[i][0:17] == 'Vertex number > 9':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            largePDBfrag.append(s)
        
largePDB = list(set(largePDB))
largePDB.remove('\n')
largePDBfrag = list(set(largePDBfrag))
largePDBfrag.remove('\n')
DualGraph_lib = {}


# Record all graphs appeared
for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            fname = id+'_'
            for x in c:
                fname += x
            # if fname not in largePDBfrag:
            graph_Collector(id, c, DualGraph_lib)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    # if fname not in largePDBfrag:
                    graph_Collector(id, c, DualGraph_lib)
                    i += 4

    
for id in largePDB:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            if i+1 < len(lines):
                c = lines[i+1].split()[1:]
                graph_Collector_largePDB(id, c, DualGraph_lib)
                i += 5
                
                while i < len(lines):
                    if lines[i] != '\n':
                        c = lines[i].split()[1:]
                        graph_Collector_largePDB(id, c, DualGraph_lib)
                        i += 4


graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))

with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library_Sub.txt', 'w') as f:
    f.write('There are in total ' +str(len(graphlist)) + ' existing graphs.\n')
    for graph in graphlist:
        DualGraph_lib[graph] = list(set(DualGraph_lib[graph]))
        hasDNA = 0
        multiChains = 0
        for fname in DualGraph_lib[graph]:
            if '*' in fname:
                hasDNA += 1
            elif '#' in fname:
                multiChains += 1
        f.write('%-15s %10s %10s# %10s*\n' % (graph+':', str(len(DualGraph_lib[graph])-hasDNA-multiChains), str(multiChains), str(hasDNA)))

                
Dual_Count = {}

with open('Dual_Library_Sub.txt', 'r') as f:
    f.readline()
    lines = f.readlines()
        
for l in lines:
    graph = l.split(':')[0]
    count = int(l.split()[1]) + int(l.split()[2].split('#')[0]) + int(l.split()[3].split('*')[0])
    if count in Dual_Count:
        Dual_Count[count] = Dual_Count[count] + graph + ', '
    else:
        Dual_Count[count] = graph + ', '
        
Dual_Count = dict(sorted(Dual_Count.items()))


with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library_Sub_Graphs.txt', 'w') as f:
    for graph in graphlist:
        f.write(graph+': ')
        for x in DualGraph_lib[graph]:
            f.write(x+', ')
        f.write('\n')

                
with open('Dual_Library_Sub_Count.txt', 'w') as f:        
    for count in Dual_Count:
        f.write(str(count) + ':\t')
        f.write(Dual_Count[count]+'\n')


    

##################### Get titles for popular graphs ############################################
PG = ['7_52', '7_814', '8_19']
PG_titles = {}

with open('Dual_Library_Sub_Graphs.txt', 'r') as f:
    lines = f.readlines()

for l in lines:
    graph = l.split(':')[0]
    if graph in PG:
        ids_raw = l.split()[1:]
        titles = []
        title = ''
        for s in ids_raw:
            id = s.split('_')[0]
            if os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb'):
                with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f1:     
                    title = f1.readline().split('\t')[3].split('\n')[0]
            elif os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif'):
                with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif', 'r') as cif:
                    ciflines = cif.readlines()
                    i = 0
                    found = False
                    while i < len(ciflines) and not found:
                        if '_struct.title' in ciflines[i]:
                            if len(ciflines[i].split()) < 2:
                                title = ciflines[i+1].split('\n')[0]
                                found = True
                        i += 1
            # with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f1:     
            #     title = f1.readline().split('\t')[3].split('\n')[0]
            titles.append(title + ' (' + s + ')')            
        PG_titles[graph] = titles
        
with open('Dual_Library_Sub_Popular.txt', 'w') as f:
    for graph in PG_titles:
        f.write('%-10s %s\n' % (graph+':', PG_titles[graph][0]))        
        for i in range(1,len(PG_titles[graph])):
            x = PG_titles[graph][i]
            f.write('%10s %s\n' % ('', x))
        f.write('\n')
        

########################### Get graph library for FSE/Riboswitch structures #################################
with open('/Users/qz886/Desktop/Dual_Library_Update/FSE_rep_list.txt', 'r') as f:
    lines = f.readlines()
    
ID = []
for l in lines:
    ID.append(l.split('\n')[0])


largePDB = []
largePDBfrag = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i][0:17] == 'Too many helices:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            largePDBfrag.append(s)
            
    if lines[i][0:17] == 'Vertex number > 9':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            largePDBfrag.append(s)
        
largePDB = list(set(largePDB))
largePDB.remove('\n')
largePDBfrag = list(set(largePDBfrag))
largePDBfrag.remove('\n')

DualGraph_lib = {}
        


############################### File writing with PDB titles #####################################
# Record all graphs appeared
for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            fname = id+'_'
            for x in c:
                fname += x
            if fname not in largePDBfrag:
                graph_Collector_title(id, c, DualGraph_lib)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    if fname not in largePDBfrag:
                        graph_Collector_title(id, c, DualGraph_lib)
                    i += 4

    
    if id in largePDB:
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
            lines = f.readlines()
            
        for i in range(len(lines)):
            if lines[i] == 'Interacting substructures\n':
                if i+1 < len(lines):
                    c = lines[i+1].split()[1:]
                    graph_Collector_largePDB_title(id, c, DualGraph_lib)
                    i += 5
                    
                    while i < len(lines):
                        if lines[i] != '\n':
                            c = lines[i].split()[1:]
                            graph_Collector_largePDB_title(id, c, DualGraph_lib)
                            i += 4


with open('/Users/qz886/Desktop/Dual_Library_Update/FSE.txt', 'r') as f:
    lines = f.readlines()

RS_type = {}
for l in lines:
    id = l.split('\t')[0]
    type = l.split('\t')[1].split('\n')[0]
    RS_type[id] = type


with open('FSE_Library_Sub_GraphTitles.txt', 'w') as f:
    if '1_1' in DualGraph_lib:
        DualGraph_lib.pop('1_1')
    graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))
    
    for graph in graphlist:
        DualGraph_lib[graph] = list(set(DualGraph_lib[graph]))
        id = DualGraph_lib[graph][0].split()[0]
        f.write(graph+'\t'+RS_type[id]+' ('+id+')\n')
        for i in range(1,len(DualGraph_lib[graph])):
            id = DualGraph_lib[graph][i].split()[0]
            f.write('\t'+RS_type[id]+' ('+id+')\n')


RS = []  
RS_type = []
for l in lines:
    RS.append(l.split('\t')[0])
    RS_type.append(l.split('\t')[1].split('\n')[0])
    
    
with open('FSE_Subgraphs.txt', 'w') as f:
    if '1_1' in DualGraph_lib:
        DualGraph_lib.pop('1_1')
    graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))
    
    for i in range(len(RS)):
        id = RS[i]
        type = RS_type[i]
        G = []
        for graph in graphlist:
            if id in DualGraph_lib[graph]:
                G.append(graph)
        G.sort()
        f.write(id+'\t'+type+'\t'+G[0]+'\t')        
        for i in range(1,len(G)):
            f.write(G[i]+'\t')
        f.write('\n')




##################################### Plot Heatmap ###################################
ribotypes = {'TPP':0, 'AqCbl':1, 'SAM-I':2, 'SAM-IV':3, 'SAM-I/IV':4, 'SAM-II':5, 'SAM-V':6, 'SAM-III':7, 'SAM-SAH':8, 'SAH':9, 'FMN':10, 'THF':11, 'ZTP':12, 'C-di-GMP-I':13, 'C-di-GMP-II':14, 'C-AMP-GMP':15, 'C-di-AMP':16, 'Glutamine':17, 'Glycine':18, 'Lysine':19, 'PreQ1-I':20, 'PreQ1-II':21, 'preQ1-III':22, 'Guanine':23, 'Adenine':24, "2'-dG-I":25, "2'-dG-II":26, 'Fluoride':27, 'Mn2+':28, 'Mg2+-I':29, 'NiCo':30, 'Guanidine-I':31, 'Guanidine-II':32, 'Guanidine-III':33, 'glmS':34}


if '1_1' in DualGraph_lib:
    DualGraph_lib.pop('1_1')
graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))

ribomatrix = np.zeros((35,len(graphlist)))

for i in range(len(graphlist)):
    graph = graphlist[i]
    DualGraph_lib[graph] = list(set(DualGraph_lib[graph]))
    for riboid in DualGraph_lib[graph]:
        id = riboid.split()[0]
        index = ribotypes[RS_type[id]]
        if index < 12: # coenzymes
            ribomatrix[index, i] = 1
        elif index < 17: # signaling molecules
            ribomatrix[index, i] = 2
        elif index < 20: # amino acids
            ribomatrix[index, i] = 3
        elif index < 27: # nucleotide
            ribomatrix[index, i] = 4
        elif index < 31: # ions
            ribomatrix[index, i] = 5
        else: # other metabolites
            ribomatrix[index, i] = 6
            
ribodf = pd.DataFrame(ribomatrix, columns=graphlist)
ribodf.index = ['TPP', 'AqCbl', 'SAM-I', 'SAM-IV', 'SAM-I/IV', 'SAM-II', 'SAM-V', 'SAM-III', 'SAM-SAH', 'SAH', 'FMN', 'THF', 'ZTP', 'C-di-GMP-I', 'C-di-GMP-II', 'C-AMP-GMP', 'C-di-AMP', 'Glutamine', 'Glycine', 'Lysine', 'PreQ1-I', 'PreQ1-II', 'preQ1-III', 'Guanine', 'Adenine', "2'-dG-I", "2'-dG-II", 'Fluoride', 'Mn2+', 'Mg2+-I', 'NiCo', 'Guanidine-I', 'Guanidine-II', 'Guanidine-III', 'glmS']

ribocolor = colors.ListedColormap(['white', 'blue', 'green', 'orange', 'red', 'purple', 'cyan'])
sns.set(rc = {'figure.figsize':(20,10)})
riboheat = sns.heatmap(ribodf, xticklabels=True, cmap = ribocolor)
figure = riboheat.get_figure()    
figure.savefig('Riboswitch_heatmap.eps')



##################################### Count substructure compositions ###################################
SingleRNA = []
MultiRNA = []
DNA = []

with open('/Users/qz886/Desktop/Dual_Library_Update/list_file.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')

# These two pdb files have unusual structures, no vertex
manual = ['7BPG', '7BPF']
for ma in manual:
    ID.remove(ma)

smallPDB = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i][0:13] == 'No base pair:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            smallPDB.append(s)
            
    if lines[i][0:10] == 'No vertex:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            smallPDB.append(s)     
        
smallPDB = list(set(smallPDB))
smallPDB.remove('\n')


for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            fname = id+'_'
            for x in c:
                fname += x
            if fname not in smallPDB:
                with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'r') as f:
                    results = f.readlines()
                if results[0] != '1_1\n':
                    if isDNA(id, c):
                        DNA.append(fname)
                    else:
                        if len(c)>1:
                            MultiRNA.append(fname)
                        else:
                            SingleRNA.append(fname)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    if fname not in smallPDB:
                        if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt'):
                            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'r') as f:
                                results = f.readlines()
                            if results[0] != '1_1\n':
                                if isDNA(id, c):
                                    DNA.append(fname)
                                else:
                                    if len(c)>1:
                                        MultiRNA.append(fname)
                                    else:
                                        SingleRNA.append(fname)
                    i += 4
                    

SingleRNA = list(set(SingleRNA))
MultiRNA = list(set(MultiRNA))
DNA = list(set(DNA))


print('Single RNA Chain: %d' %len(SingleRNA))
print('Multiple RNA Chains: %d' %len(MultiRNA))
print('DNA Chains: %d' %len(DNA))
print('Total: %d' %(len(SingleRNA)+len(MultiRNA)+len(DNA)))