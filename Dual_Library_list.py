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



def graph_Collector(id, chains, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs):
    fname = id+'_'
    chainList = []
    for x in chains:
        fname += x
        chainList.append(x.split('-')[0])
    chainList = list(set(chainList))
                
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
            
        if lines[0] == 'No base pair':
            NoBP.append(fname)
            
        elif lines[0] == 'Too many helices!':
            LargeHelices.append(fname)
            pass
            
        elif lines[0] == 'Vertex number > 9\n':
            LargeGraph.append(fname)
            pass
            
        elif lines[0] == 'No vertex\n':
            NoVertex.append(fname)
            
        elif lines[0] == 'No matching graph exists (even if the vertex number is between 2 and 9).':
            NoGraph.append(fname)
        
        elif '_' in lines[0]: # this RNA gets assigned of one dual graph
            AllGraphs.append(fname)
            graph = lines[0].split()[0]
            if graph in DualGraph_lib:
                if isDNA(id, chains):
                    DualGraph_lib[graph].append(fname+'*')
                else:
                    if len(chainList)>1:
                        DualGraph_lib[graph].append(fname+'#')
                    else:
                        DualGraph_lib[graph].append(fname)
            else:
                if isDNA(id, chains):
                    DualGraph_lib[graph] = [fname+'*']
                else:
                    if len(chainList)>1:
                        DualGraph_lib[graph] = [fname+'#']
                    else:
                        DualGraph_lib[graph] = [fname]

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
            
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        if '_' in lines[0]: # this RNA gets assigned of one dual graph
            graph = lines[0].split()[0]
            if graph in DualGraph_lib:
                DualGraph_lib[graph].append(title+' ('+fname+')')
            else:
                DualGraph_lib[graph] = [title+' ('+fname+')']
    else:
        print('Cannot identify '+fname+'.txt')



def graph_Collector_largePDB(id, chains, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs):
    fname = id+'_'
    chainList = []
    for x in chains:
        fname += x
        chainList.append(x.split('-')[0])
    chainList = list(set(chainList))
                
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
            
        if lines[0] == 'No base pair':
            NoBP.append(fname)
            
        elif lines[0] == 'Too many helices!':
            LargeHelices.append(fname)
            
        elif lines[0] == 'Vertex number > 9\n':
            LargeGraph.append(fname)
            
        elif lines[0] == 'No vertex\n':
            NoVertex.append(fname)
            
        elif lines[0] == 'No matching graph exists (even if the vertex number is between 2 and 9).':
            NoGraph.append(fname)
        
        elif '_' in lines[0]: # this RNA gets assigned of one dual graph
            AllGraphs.append(fname)
            graph = lines[0].split()[0]
            if graph in DualGraph_lib:
                if isDNA_largePDB(id, chains):
                    DualGraph_lib[graph].append(fname+'*')
                else:
                    if len(chainList)>1:
                        DualGraph_lib[graph].append(fname+'#')
                    else:
                        DualGraph_lib[graph].append(fname)
            else:
                if isDNA_largePDB(id, chains):
                    DualGraph_lib[graph] = [fname+'*']
                else:
                    if len(chainList)>1:
                        DualGraph_lib[graph] = [fname+'#']
                    else:
                        DualGraph_lib[graph] = [fname]

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
            
    if os.path.isfile('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt'):
        with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt', 'r') as f:
            lines = f.readlines()
        
        if '_' in lines[0]: # this RNA gets assigned of one dual graph
            graph = lines[0].split()[0]
            if graph in DualGraph_lib:
                DualGraph_lib[graph].append(title+' ('+fname+')')
            else:
                DualGraph_lib[graph] = [title+' ('+fname+')']
    else:
        print('Cannot identify '+fname+'.txt')
        



########################## Get graph library for all structures #################################
with open('/Users/qz886/Desktop/Dual_Library_Update/list_file.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')

# These two pdb files have unusual structures, no vertex
manual = ['7BPG', '7BPF']
for ma in manual:
    ID.remove(ma)

largePDB = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i][0:17] == 'Too many helices:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            
    if lines[i][0:17] == 'Vertex number > 9':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
        
largePDB = list(set(largePDB))
largePDB.remove('\n')

        
NoBP = []
LargeHelices = []
LargeGraph = []
NoVertex = []
NoGraph = []
AllGraphs = []
DualGraph_lib = {}


# Record all graphs appeared
for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            graph_Collector(id, c, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    graph_Collector(id, c, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs)
                    i += 4
    
    if id in largePDB:
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
            lines = f.readlines()
            
        for i in range(len(lines)):
            if lines[i] == 'Interacting substructures\n':
                if i+1 < len(lines):
                    c = lines[i+1].split()[1:]
                    graph_Collector_largePDB(id, c, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs)
                    i += 5
                    
                    while i < len(lines):
                        if lines[i] != '\n':
                            c = lines[i].split()[1:]
                            graph_Collector_largePDB(id, c, NoBP, LargeHelices, LargeGraph, NoVertex, NoGraph, DualGraph_lib, AllGraphs)
                            i += 4

        
with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'w') as f:
    f.write('No base pair: '+ str(len(NoBP)) +'\n')
    for id in NoBP:
        f.write(id+',\t')
    f.write('\n\n\n')
    
    f.write('Too many helices: '+ str(len(LargeHelices)) +'\n')
    for id in LargeHelices:
        f.write(id+',\t')
    f.write('\n\n\n')
        
    f.write('Vertex number > 9: '+ str(len(LargeGraph)) +'\n')
    for id in LargeGraph:
        f.write(id+',\t')
    f.write('\n\n\n')
        
    f.write('No vertex: '+ str(len(NoVertex)) +'\n')
    for id in NoVertex:
        f.write(id+',\t')
    f.write('\n\n\n')
        
    f.write('No matching graph exists (even if the vertex number is between 2 and 9): '+ str(len(NoGraph)) +'\n')
    for id in NoGraph:
        f.write(id+',\t')
    f.write('\n\n\n')

# DualGraph_lib.pop('1_1')
graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library.txt', 'w') as f:
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

        
with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_Graphs.txt', 'w') as f:
    for graph in graphlist:
        f.write(graph+': ')
        for x in DualGraph_lib[graph]:
            f.write(x+', ')
        f.write('\n')



Dual_Count = {}

with open('Dual_Library.txt', 'r') as f:
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
                
with open('Dual_Library_Count.txt', 'w') as f:        
    for count in Dual_Count:
        f.write(str(count) + ':\t')
        f.write(Dual_Count[count]+'\n')



###################### Get titles for popular graphs ############################################
PG = ['5_2']
PG_titles = {}

with open('Dual_Library_Graphs.txt', 'r') as f:
    lines = f.readlines()

for l in lines:
    graph = l.split(':')[0]
    if graph in PG:
        ids_raw = l.split()[1:]
        titles = []
        for s in ids_raw:
            id = s.split('_')[0]
            with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f1:     
                title = f1.readline().split('\t')[3].split('\n')[0]
            titles.append(title + ' (' + s + ')')            
        PG_titles[graph] = titles
        
with open('Dual_Library_Popular.txt', 'w') as f:
    for graph in PG_titles:
        f.write('%-10s %s\n' % (graph+':', PG_titles[graph][0]))        
        for i in range(1,len(PG_titles[graph])):
            x = PG_titles[graph][i]
            f.write('%10s %s\n' % ('', x))
        f.write('\n')
        


###################### Find organisms for 7_1311 and 6_263, 5S rRNA ##############################     
PG = ['4_27']
PG_titles = {}

with open('Dual_Library_Graphs.txt', 'r') as f:
    lines = f.readlines()

for l in lines:
    graph = l.split(':')[0]
    if graph in PG:
        ids_raw = l.split()[1:]
        titles = []
        for s in ids_raw:
            id = s.split('_')[0]
            org = ''
            # see if PDB or CIF file
            if os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif'):
                with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif', 'r') as cif:
                    ciflines = cif.readlines()
                    for ind in range(len(ciflines)):
                        if ciflines[ind] == '_entity_src_nat.details \n':
                            orgname = re.findall(r'\d+ ([?a-zA-Z \']+\'[ -?\w\']+) \d+', ciflines[ind+1])
                            if len(orgname) > 0:
                                org = orgname[0]
            elif os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb'):
                with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb', 'r') as cif:
                    ciflines = cif.readlines()
                    for cl in ciflines:
                        if 'ORGANISM_COMMON' in cl:
                            org = cl.split('ORGANISM_COMMON:')[1].split('\n')[0]
            titles.append(org + ' (' + id + ')')            
        PG_titles[graph] = titles
        
with open('5S rRNA.txt', 'w') as f:
    for graph in PG_titles:
        f.write('%-10s %s\n' % (graph+':', PG_titles[graph][0]))        
        for i in range(1,len(PG_titles[graph])):
            x = PG_titles[graph][i]
            f.write('%10s %s\n' % ('', x))
        f.write('\n')



########################### Get graph library for all Riboswitch or FSE (change name) RNA structures #################################
with open('/Users/qz886/Desktop/Dual_Library_Update/Riboswitch_rep_list.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')


with open('/Users/qz886/Desktop/Dual_Library_Update/Riboswitch.csv', 'r') as f:
    lines = f.readlines()

RS_type = {}
for l in lines:
    id = l.split(',')[0]
    type = l.split(',')[1]
    RS_type[id] = type
    


largePDB = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library_noGraph.txt', 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i][0:17] == 'Too many helices:':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
            
    if lines[i][0:17] == 'Vertex number > 9':
        IDstring = lines[i+1].split(',\t')
        for s in IDstring:
            largePDB.append(s.split('_')[0])
       
largePDB = list(set(largePDB))
largePDB.remove('\n')

        
NoBP = []
LargeHelices = []
LargeGraph = []
NoVertex = []
NoGraph = []
DualGraph_lib = {}



# Record all graphs appeared
for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            graph_Collector_title(id, c, DualGraph_lib)
            i += 5
           
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
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


with open('Riboswitch_graphs.txt', 'w') as f:
    if '1_1' in DualGraph_lib:
        DualGraph_lib.pop('1_1')
    graphlist = sorted(DualGraph_lib, key=lambda s: list(map(int, s.split('_'))))
   
    for graph in graphlist:
        DualGraph_lib[graph] = list(set(DualGraph_lib[graph]))
        
        id = DualGraph_lib[graph][0].split('(')[1].split(')')[0].split('_')[0]
        f.write('%-10s %s\n' % (graph+':', RS_type[id]+' ('+id+')'))        
        for i in range(1,len(DualGraph_lib[graph])):
            id = DualGraph_lib[graph][i].split('(')[1].split(')')[0].split('_')[0]
            f.write('%10s %s\n' % ('', RS_type[id]+' ('+id+')'))
        f.write('\n')
