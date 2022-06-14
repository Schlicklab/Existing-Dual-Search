#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 22:30:00 2022

@author: qiyaozhu
"""

import os


with open('Dual_Library_Old.txt', 'r') as f:
    f.readline()
    old = f.readlines()
    
with open('Dual_Library.txt', 'r') as f:
    f.readline()
    new = f.readlines()
    
G_old = []
G_new = []

for l in old:
    G_old.append(l.split(':')[0])
    
for l in new:
    G_new.append(l.split(':')[0])

Old_only = []
New_only = []

for i in range(len(old)):
    if G_old[i] not in G_new:
        Old_only.append(old[i])

for i in range(len(new)):
    if G_new[i] not in G_old:
        New_only.append(new[i])


with open('Dual_Library_Graphs.txt', 'r') as f:
    graphlist = f.readlines()

GraphList = {}
for gl in graphlist:
    GraphList[gl.split(':')[0]] = gl.split(':')[1]

        
with open('Dual_Library_Compare.txt', 'w') as f:
    f.write('Old Library only '+str(len(Old_only))+' graphs:\n')
    for l in Old_only:
        f.write(l)
        
    f.write('\n\nNew Library only '+str(len(New_only))+' graphs:\n')
    for l in New_only:
        f.write(l.split('\n')[0]+'\t\t')
        f.write(GraphList[l.split(':')[0]])


New_only_Count = {}
for l in New_only:
    graph = l.split('\n')[0].split(':')[0]
    count = int(l.split('\n')[0].split()[1]) + int(l.split('\n')[0].split()[2].split('#')[0]) + int(l.split('\n')[0].split()[3].split('*')[0])
    if count in New_only_Count:
        New_only_Count[count] = New_only_Count[count] + graph + ', '
    else:
        New_only_Count[count] = graph + ', '
        
New_only_Count = dict(sorted(New_only_Count.items()))
                
with open('Dual_Library_Compare_Count.txt', 'w') as f:        
    f.write('New Library only '+str(len(New_only))+' graphs:\n')
    for count in New_only_Count:
        f.write(str(count) + ':\t')
        f.write(New_only_Count[count]+'\n')
    

