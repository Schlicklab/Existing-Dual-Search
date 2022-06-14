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


def readDualGraphs():

    DualGraphsLib=[]
    Graphs=[]
    DualGraphsLib.append(Graphs) # for vertex == 1, no graphs are there
    for i in range(2,10): # will read dual graphs from 2-9 vertices (10 used as the range function stops before the last number)

        Graphs=[]
        file_eigen = "/Users/qz886/Desktop/Dual_Graph_Codes/DualEig/%dEigen"%i
        file_adjMat = "/Users/qz886/Desktop/Dual_Graph_Codes/DualAdj/V%dAdjDG"%i
        
        loadEigenvalues(Graphs,i,file_eigen) # load eigen values for dual graphs for vertex number i
        loadAdjMatrices(Graphs,i,file_adjMat) # load adjacency matrices for dual graphs for vertex number i

        DualGraphsLib.append(Graphs)
    
    return DualGraphsLib

              

def get_Subgraphs(A, DualGraphsLib):
    
    with open('Adj.txt', 'w') as f:
        n = len(A)
        for i in range(0,n):    
            for j in range(0,n):   
                f.write(str(A[i][j])+'\t')
            f.write('\n')
    
    os.system('./dualgraph.out -input Adj.txt -len '+str(n)+' -output Sub.txt -all Sub_all.txt')
    
    with open('Sub_all.txt', 'r') as f:
        lines = f.readlines()
        
    os.system('rm -rf Adj.txt Sub.txt Sub_all.txt')
    
    subgraphs = []
    smallfails = 0
    largefails = 0
    
    for l in lines:
        if l[0] == '(':
            edges = [x.strip() for x in l.split('-')]
            numbers = [int(x) for x in re.findall(r"[\w']+",l)]
            vertices =list(set(numbers))
            
            matrix = []
            for a in range(0,len(vertices)):
                tempArray = []
                for b in range(0,len(vertices)):
                    tempArray.append(0)
                matrix.append(tempArray)
            
            for j in range(len(edges)-1): #last entry is empty as there is a "-" at the end of each line
                edge = edges[j]
                indices = [int(x) for x in re.findall(r"[\w']+",edge)] #read each pair (edge). i.e.  (11,3)
                m = vertices.index(indices[0]) #determine the order (index) of the first vertex of the edge. For (11,3) it is 11 and the index is 3 according to vertices
                n = vertices.index(indices[1]) ##determine the order (index) of the first vertex of the edge. It is 3 and the index is 0               
                matrix[m][n]+=1 #increase the number of connections in the adjacency matrix. matrix[3][0] will be increased 1
                matrix[n][m]+=1 #since the matrix is symmetric, increase matrix[0][3] 1.
                
            N=len(vertices)
            if N==1:
                print("1_1\n")
                smallfails += 1
            elif N>9:
                print('Vertex number > 9\n')
                largefails += 1
            else:
                eigen = calcEigenValues(matrix) # calculate the eigen values for the subgraph matrix
                subgraphID = searchtoAssignID(DualGraphsLib[N-1],0,len(DualGraphsLib[N-1])-1,eigen,matrix)
                subgraphs.append(subgraphID)
                print(subgraphID)
    
    return subgraphs, smallfails, largefails



def ctToSubgraphs(fname, DualGraphsLib):
    if os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct"):
        RNA = getCTInfo("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct") 
        bpexist = countHelices(RNA)
        
        if bpexist:
            print('Extracting dual graphs for '+fname)
            changeHelices(RNA)
            RNA.makeMatrices()
            connectHelices(RNA)
            vertexOrder = []
            for i in range(0,len(RNA.adjMatrix)):
                vertexOrder.append(0)
            correctHNumbers(RNA)
                
#            DualGraphsLib = readDualGraphs()    
            subgraphs, smallfails, largefails = get_Subgraphs(RNA.adjMatrix, DualGraphsLib)
            subgraphs = sorted(list(set(subgraphs)))
                
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt', 'w') as f:
                if subgraphs: # not empty
                    for g in subgraphs:
                        f.write(g+'\n')
                else:
                    f.write('No subgraph for 1_1\n')
                if smallfails:
                    f.write('Vertex number < 2\n')
                if largefails:
                    f.write('Vertex number > 9\n')        
        else:
            print(fname+' has no base pair!')
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub/'+fname+'.txt', 'w') as f:
                f.write('No base pair')

                
def ctToSubgraphs_largePDB(fname, DualGraphsLib):
    if os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct"):
        RNA = getCTInfo("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct") 
        bpexist = countHelices(RNA)
        
        if bpexist:
            print('Extracting dual graphs for '+fname)
            changeHelices(RNA)
            RNA.makeMatrices()
            connectHelices(RNA)
            vertexOrder = []
            for i in range(0,len(RNA.adjMatrix)):
                vertexOrder.append(0)
            correctHNumbers(RNA)
                
#            DualGraphsLib = readDualGraphs()    
            subgraphs, smallfails, largefails = get_Subgraphs(RNA.adjMatrix, DualGraphsLib)
            subgraphs = sorted(list(set(subgraphs)))
                
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt', 'w') as f:
                if subgraphs: # not empty
                    for g in subgraphs:
                        f.write(g+'\n')
                else:
                    f.write('No subgraph for 1_1\n')
                if smallfails:
                    f.write('Vertex number < 2\n')
                if largefails:
                    f.write('Vertex number > 9\n')        
        else:
            print(fname+' has no base pair!')
            with open('/Users/qz886/Desktop/Dual_Graph_Codes/PDB_DSSR_Dual_Sub_largePDB/'+fname+'.txt', 'w') as f:
                f.write('No base pair')




with open('/Users/qz886/Desktop/Dual_Library_Update/list_file.txt', 'r') as f:
    lines = f.readlines()
ID = lines[0].split(',')

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


DualGraphsLib=[]
DualGraphsLib = readDualGraphs()


## Get all subgraphs
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
                ctToSubgraphs(fname, DualGraphsLib)
            i += 5
            
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    if fname not in largePDBfrag:
                        ctToSubgraphs(fname, DualGraphsLib)
                    i += 4
        
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
                ctToSubgraphs_largePDB(fname, DualGraphsLib)
                i += 5
                
                while i < len(lines):
                    if lines[i] != '\n':
                        c = lines[i].split()[1:]
                        fname = id+'_'
                        for x in c:
                            fname += x
                        ctToSubgraphs_largePDB(fname, DualGraphsLib)
                        i += 4
