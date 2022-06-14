#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 19:59:28 2021

@author: qz886
"""

from igraph import *
from ClassesFunctions import *
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
        file_eigen = "DualEig/%dEigen"%i
        file_adjMat = "DualAdj/V%dAdjDG"%i
        
        loadEigenvalues(Graphs,i,file_eigen) # load eigen values for dual graphs for vertex number i
        loadAdjMatrices(Graphs,i,file_adjMat) # load adjacency matrices for dual graphs for vertex number i

        DualGraphsLib.append(Graphs)
    
    return DualGraphsLib
    

def get_Adj(ID):
    
    n = int(ID.split('_')[0])
    Graphs = []
    eigen_file = "DualEig/%dEigen"%n
    adj_file = "DualAdj/V%dAdjDG"%n
    loadEigenvalues(Graphs,n,eigen_file)
    loadAdjMatrices(Graphs,n,adj_file)
    
    for g in Graphs:
        if g.graphID == ID:
            A = g.adjMatrix
            return A
    

def plotGraph(ID):
    
    n = int(ID.split('_')[0])
    A = get_Adj(ID)
    g = Graph()
    g.add_vertices(n)

    mylayout=g.layout_circle()
    bbox=BoundingBox(400,400)

    for i in range(0,n):    
        for j in range(i,n):   
            for k in range(0,int(A[i][j])):
                g.add_edge(i,j)
    if n > 3:
        mylayout=g.layout_fruchterman_reingold()
    
    plot(g, "Subgraphs/%s.png"%ID,vertex_color="#b5651d", vertex_frame_color="#b5651d", vertex_size=40, edge_color="#b5651d", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
                

def get_Subgraphs(ID):
    
    n = int(ID.split('_')[0])
    A = get_Adj(ID)
    with open(ID+'Adj.txt', 'w') as f:
        for i in range(0,n):    
            for j in range(0,n):   
                f.write(str(A[i][j])+'\t')
            f.write('\n')
    
    os.system('./dualgraph.out -input '+ID+'Adj.txt -len '+str(n)+' -output '+ID+'Sub.txt -all '+ID+'Sub_all.txt')
    
    with open(ID+'Sub_all.txt', 'r') as f:
        lines = f.readlines()
        
    os.system('rm -rf '+ID+'Adj.txt '+ID+'Sub.txt '+ID+'Sub_all.txt')
    
    DualGraphsLib = readDualGraphs()
    subgraphs = []
    
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
            elif N>9:
                print('Vertex number > 9\n')
            else:
                eigen = calcEigenValues(matrix) # calculate the eigen values for the subgraph matrix
                subgraphID = searchtoAssignID(DualGraphsLib[N-1],0,len(DualGraphsLib[N-1])-1,eigen,matrix)
                subgraphs.append(subgraphID)
                print(subgraphID)
    
    return subgraphs


def is_pseudoknot(ID):
    
    n = int(ID.split('_')[0])
    A = get_Adj(ID)
    with open(ID+'Adj.txt', 'w') as f:
        for i in range(0,n):    
            for j in range(0,n):   
                f.write(str(A[i][j])+'\t')
            f.write('\n')
    
    os.system('./dualgraph.out -input '+ID+'Adj.txt -len '+str(n)+' -output '+ID+'Sub.txt -all '+ID+'Sub_all.txt')
    
    with open(ID+'Sub.txt', 'r') as f:
        lines = f.readlines()
    
    is_pknot = None
    
    for l in lines:
        if 'number of PK blocks:' in l:
            pkn = int(l.split('\n')[0].split(': ')[1])
            if pkn > 0:
                is_pknot = True
            elif pkn == 0:
                is_pknot = False
            else:
                print('Error in finding pseudoknots!')
            break
        
    os.system('rm -rf '+ID+'Adj.txt '+ID+'Sub.txt '+ID+'Sub_all.txt')
    
    return is_pknot

                
# ID = sys.argv[1]
# subgraphs = get_Subgraphs(ID)

# for id in subgraphs:
#     plotGraph(id)
    
    
    
    
    