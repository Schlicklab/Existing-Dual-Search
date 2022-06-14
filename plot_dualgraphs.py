from igraph import *
from ClassesFunctions import *
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as patches
from Subgraphs import *


def plotDual(g):
    
    A = np.array(g.adjMatrix)
    g1 = Graph()
    n = A[0].size
    g1.add_vertices(n)
        
    mylayout=g1.layout_circle()
    bbox=BoundingBox(400,400)

    for i in range(0,n):    
        for j in range(i,n):   
            for k in range(0,int(A[i][j])):
                g1.add_edge(i,j)
    
    mylayout=g1.layout_fruchterman_reingold()
                
    # vlabel = np.argsort(FV)+1
    vlsize = 30
    vldist = 1.5
    
    plot(g1, "%s.eps"%g.graphID,vertex_color="#FF0000", vertex_frame_color="#FF0000", vertex_size=40, edge_color="#FF0000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
 

n = 5

Graphs = []
eigen_file = "DualEig/%dEigen"%n
adj_file = "DualAdj/V%dAdjDG"%n
loadEigenvalues(Graphs,n,eigen_file)
loadAdjMatrices(Graphs,n,adj_file)

for g in Graphs:
    print(g.graphID)
    if g.graphID in ['5_96']:
        plotDual(g)





###################### Plot the composite figure of all existing dual graphs, with markers ######################

def plotGraph(IDs, path):
    
    columns = 14
    rows = int(ceil(len(IDs)/columns))
    fig = plt.figure(figsize=(20*columns, 25*rows))
    
    Vn = []
    for id in IDs:
        Vn.append(int(id.split('_')[0]))
    Vn = sorted(list(set(Vn)))
    
    AllGraphs = {}
    for n in Vn:
        Graphs = []
        eigen_file = "DualEig/%dEigen"%n
        adj_file = "DualAdj/V%dAdjDG"%n
        loadEigenvalues(Graphs,n,eigen_file)
        loadAdjMatrices(Graphs,n,adj_file)
        AllGraphs[n] = Graphs
    
    for id in IDs:
        n = int(id.split('_')[0])
        Graphs = AllGraphs[n]
        for g in Graphs:
            if g.graphID == id:
                A = g.adjMatrix
                g1 = Graph()
                g1.add_vertices(n)
        
                mylayout=g1.layout_circle()
                bbox=BoundingBox(400,400)

                for i in range(0,n):    
                    for j in range(i,n):   
                        for k in range(0,int(A[i][j])):
                            g1.add_edge(i,j)
                if n > 3:
                    mylayout=g1.layout_fruchterman_reingold()
                
                isOld = IDs[id][1]
                mark = IDs[id][2]
                if mark == 'R':
                    if isOld:
                        plot(g1, path+id+'.png',vertex_color="#FF0000", vertex_frame_color="#FF0000", vertex_size=40, edge_color="#FF0000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
                    else:
                        plot(g1, path+id+'.png',vertex_color="#FF0000", vertex_frame_color="#FF0000", vertex_size=40, edge_color="#FF0000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="#FFFF00", margin=(75,75,75,75))
                elif mark == 'M':
                    if isOld:
                        plot(g1, path+id+'.png',vertex_color="#BF40BF", vertex_frame_color="#BF40BF", vertex_size=40, edge_color="#BF40BF", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
                    else:
                        plot(g1, path+id+'.png',vertex_color="#BF40BF", vertex_frame_color="#BF40BF", vertex_size=40, edge_color="#BF40BF", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="#FFFF00", margin=(75,75,75,75))
                elif mark == 'D':
                    if isOld:
                        plot(g1, path+id+'.png',vertex_color="#008000", vertex_frame_color="#008000", vertex_size=40, edge_color="#008000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="white", margin=(75,75,75,75))
                    else:
                        plot(g1, path+id+'.png',vertex_color="#008000", vertex_frame_color="#008000", vertex_size=40, edge_color="#008000", edge_width=4, layout=mylayout, rescale=True, bbox=bbox, background="#FFFF00", margin=(75,75,75,75))
                else:
                    print('MARKER READING PROBLEM!!!!')
    
    pos = 1
    for id in IDs:
        img = mpimg.imread(path+id+'.png')
        subf = fig.add_subplot(rows, columns, pos)
        is_pknot = is_pseudoknot(id)
        round = IDs[id][0]
        if round == 1:
            subf.set_title(id, size=180)
        else:
            subf.set_title(id+' *', size=180)
        if is_pknot:
            rect = patches.Rectangle((15, 15), 370, 370, linewidth=12, edgecolor='k', facecolor='none')
            subf.add_patch(rect)
        pos += 1
        plt.imshow(img)   
        plt.axis('off')
    plt.savefig(path+"ExistingDuals.eps")



# dual graphs found in the the first search round
IDs_first = []
with open('/Users/qz886/Desktop/Dual_Graph_Codes/firstRound/Dual_Library.txt', 'r') as f:
    f.readline()
    lines = f.readlines()
for l in lines:
    IDs_first.append(l.split(':')[0])

# dual graphs found in the the prior study   
IDs_old = []
with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library_Old.txt', 'r') as f:
    f.readline()
    old = f.readlines()
for o in old:
    IDs_old.append(o.split(':')[0])

  
IDs = {}
DNA = 0
multiRNA = 0
singleRNA = 0
new = 0

with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library.txt', 'r') as f:
    f.readline()
    lines = f.readlines()
    
for l in lines:
    id = l.split(':')[0]
    
    if id in IDs_first:
        round = 1
    else:
        round = 2
    
    if id in IDs_old:
        isOld = True
    else:
        isOld = False
        new += 1
        
    if int(l.split()[1]) > 0:
        mark = 'R'
        singleRNA += 1
    elif int(l.split()[2].split('#')[0]) > 0:
        mark = 'M'
        multiRNA += 1
    elif int(l.split()[3].split('*')[0]) > 0:
        mark = 'D'
        DNA += 1
    else:
        print('MARKER WRITING PROBLEM!!!!')
    IDs[id] = [round, isOld, mark]


path = '/Users/qz886/Desktop/Dual_Graph_Codes/Figures/'
plotGraph(IDs, path)




