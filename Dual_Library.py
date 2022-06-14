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
from plot_dualgraphs import *


def get_Graphs(fname, DualGraph_lib, NoCT):

    if os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct"):
        RNA = getCTInfo("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct")
        bpexist = countHelices(RNA)
        
        if bpexist:
            if len(RNA.Helices)>100:
                print('Too many helices!')
                with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'w') as f:
                    f.write('Too many helices!')
            
            else:                
                print('Extracting dual graphs for '+fname)
                changeHelices(RNA)
                RNA.makeMatrices()
                connectHelices(RNA)
                
                print ("Number of Vertices: " + str(len(RNA.Helices)-1))
                RNA.printAdj()
                RNA.printDeg()
                
                vertexOrder = []
                for i in range(0,len(RNA.adjMatrix)): # S.J. 07/11/2018 - to keep track of vertexOrder
                    vertexOrder.append(0)
                
                success, graph=calcEigen(RNA, vertexOrder)
                correctHNumbers(RNA)
                    
                if success == 1:
                    if graph in DualGraph_lib:
                        DualGraph_lib[graph].append(fname)
                    else:
                        DualGraph_lib[graph] = [fname]
                    
                with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'w') as f:
                    if success == 1:
                        f.write(graph+'\n')
                    elif len(RNA.adjMatrix)==1:
                        f.write('Vertex number < 2\n')
                    elif len(RNA.adjMatrix)>9:
                        f.write('Vertex number > 9\n')
                    elif len(RNA.adjMatrix)==0:
                        f.write('No vertex\n')
                    else:
                        f.write("No matching graph exists (even if the vertex number is between 2 and 9).")
    
        else:
            print(fname+' has no base pair!')
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual/'+fname+'.txt', 'w') as f:
                f.write('No base pair')
    
    else:
        NoCT.append(fname)
        


def get_Graphs_largePDB(fname, DualGraph_lib, NoCT):

    if os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct"):
        RNA = getCTInfo("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct")
        bpexist = countHelices(RNA)
        
        if bpexist:
            if len(RNA.Helices)>100:
                print('Too many helices!')
                with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt', 'w') as f:
                    f.write('Too many helices!')
            
            else:                
                print('Extracting dual graphs for '+fname)
                changeHelices(RNA)
                RNA.makeMatrices()
                connectHelices(RNA)
                
                print ("Number of Vertices: " + str(len(RNA.Helices)-1))
                RNA.printAdj()
                RNA.printDeg()
                
                vertexOrder = []
                for i in range(0,len(RNA.adjMatrix)): # S.J. 07/11/2018 - to keep track of vertexOrder
                    vertexOrder.append(0)
                
                success, graph=calcEigen(RNA, vertexOrder)
                correctHNumbers(RNA)
                    
                if success == 1:
                    if graph in DualGraph_lib:
                        DualGraph_lib[graph].append(fname)
                    else:
                        DualGraph_lib[graph] = [fname]
                    
                with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt', 'w') as f:
                    if success == 1:
                        f.write(graph+'\n')
                    elif len(RNA.adjMatrix)==1:
                        f.write('Vertex number < 2\n')
                    elif len(RNA.adjMatrix)>9:
                        f.write('Vertex number > 9\n')
                    elif len(RNA.adjMatrix)==0:
                        f.write('No vertex\n')
                    else:
                        f.write("No matching graph exists (even if the vertex number is between 2 and 9).")
    
        else:
            print(fname+' has no base pair!')
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_largePDB/'+fname+'.txt', 'w') as f:
                f.write('No base pair')
    
    else:
        NoCT.append(fname)





##################################### First round: Just get all graphs ########################################
with open('/Users/qz886/Desktop/Dual_Library_Update/MissedID_list.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')

# # These two pdb files have unusual structures, no vertex
# manual = ['7BPG', '7BPF']
# for ma in manual:
#     ID.remove(ma)

DualGraph_lib = {}        
NoCT = []        

for id in ID:
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            fname = id+'_'
            for x in c:
                fname += x
            get_Graphs(fname, DualGraph_lib, NoCT)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    get_Graphs(fname, DualGraph_lib, NoCT)
                    i += 4
    
print(NoCT)






################################## Second round: Redo graphs for large structures that have vertices >9 #######################
largePDB = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library_noGraph.txt', 'r') as f:
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

print(len(largePDB))

DualGraph_lib = {}        
NoCT = []  

        
for id in largePDB:   
            
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
        
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            if i+1 < len(lines):
                c = lines[i+1].split()[1:]
                fname = id+'_'
                for x in c:
                    fname += x
                get_Graphs_largePDB(fname, DualGraph_lib, NoCT)
                i += 5
                
                while i < len(lines):
                    if lines[i] != '\n':
                        c = lines[i].split()[1:]
                        fname = id+'_'
                        for x in c:
                            fname += x
                        get_Graphs_largePDB(fname, DualGraph_lib, NoCT)
                        i += 4
    

print(NoCT)





