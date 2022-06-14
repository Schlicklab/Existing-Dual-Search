#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 19:15:24 2021

@author: qiyaozhu
"""

import os
import os.path
from difflib import SequenceMatcher
import copy


with open('/Users/qz886/Desktop/Dual_Library_Update/list_file.txt', 'r') as f:
    lines = f.readlines()
    
ID = lines[0].split(',')
Problems = []


##################### Extract all structures from the PDB files using DSSR-3DNA ##########################

Fail = []
for id in ID:
    PDBexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/' + id + '.pdb')
    CIFexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/' + id + '.cif')
    
    if PDBexist:
        print('Extracting '+id+".pdb")
        os.system('./x3dna-dssr -i=/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb -o=/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR/'+id+'.out')
    
    elif CIFexist:
        print('Extracting '+id+".cif")
        os.system('./x3dna-dssr -i=/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif -o=/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR/'+id+'.out')
    
    else:
        print('No '+id)
        Fail.append(id)
        continue
        
print('Failed chains are:\n')
for f in Fail:
    print(f)



##################### Separate the pdb files into interacting substructures, and get 2D structure for each ##################
# These two pdb files have unusual structures, no vertex
manual = ['7BPG', '7BPF']
for ma in manual:
    ID.remove(ma)


aaList = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

Dimers = []


for id in ID:
  
    print(id)

    name = ''
    seq = ''
    db = ''
    chains = []
    chainType = {}
    subchainPos = {} # for each chain, record which subchain every residue belongs to
    chainBind = {}
    BasePairs = {} # record base pairs between every two chains
    chainSeq = {}
    chainDB = {}
    metals = {}
    protein = 'NO'
    title = ''
    pseudo = 'NO'
    NT = [] # all nucleotide information
  
  
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR/'+id+'.out', 'r') as f:
        lines = f.readlines()
  
    n = 0    
    while n<len(lines):
      
        if 'no. of metals:' in lines[n]:
            me = lines[n].split()[-1]
            if me != '0':
                me = me.split('[')[1].split(']')[0].split(',')
                for m in me:
                    metals[m.split('=')[0]] = m.split('=')[1]
            n += 1
              
        elif lines[n] == 'Secondary structures in dot-bracket notation (dbn) as a whole and per chain\n':
            name = lines[n+1]
            seq = lines[n+2]
            db = lines[n+3]
            if '[' in db:
                pseudo = 'YES'
            n += 4            
            ############### New change, now consider broken subchains as individual chains ##################
            while lines[n][0] == '>':
                c = lines[n].split('>')[1].split(' ')[0].split('-')[1]
                subchains = lines[n+1].split('\n')[0].split('&')
                subchainDB = lines[n+2].split('\n')[0].split('&')
                if len(subchains) == 1:
                    chains.append(c)
                    chainType[c] = lines[n].split(' ')[-1].split('\n')[0]
                    subchainPos[c] = [c for s in subchains[0]]
                    chainSeq[c] = subchains[0]
                    chainDB[c] = subchainDB[0]
                else:
                    subchainPos[c] = []
                    for q in range(len(subchains)):
                        chains.append(c+'-'+str(q+1))
                        chainType[c+'-'+str(q+1)] = lines[n].split(' ')[-1].split('\n')[0]
                        subchainPos[c] += [c+'-'+str(q+1) for s in subchains[q]]
                        chainSeq[c+'-'+str(q+1)] = subchains[q]
                        chainDB[c+'-'+str(q+1)] = subchainDB[q]
                n += 3
              
        elif lines[n][0:30] == 'Summary of structural features':
            n += 8
            while lines[n] != '\n':
                NT.append(lines[n])
                n += 1
              
        else:
            n +=1


    # Write the original whole structure in ct format to help later chain reorganization
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/'+id+'.out', 'w') as f:
        f.write('>'+id+'\n')
        f.write(seq.replace('&',''))
        f.write(db.replace('&',''))
    os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+id+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+id+".ct" )
  

    #################### New change, now find all (sub)chain interactions using the ct file ########################
    # reorder the chains and relabel each residue because of subchains
    reorder = []
    relabel = {}
    missingchain = ''
    respos = 1
    for nt in NT:
        cc = nt.split()[3].split('.')[0]
        ################ New change, find DSSR files have errors, some chains are not listed in the Secondary structures section, e.g. 3AVT #############################
        if cc in subchainPos:
            subcc = subchainPos[cc].pop(0)
            if subcc not in reorder:
                reorder.append(subcc)
            relabel[respos] = subcc
        else:
            if cc != missingchain:
                reorder.append(cc)
                chainSeq[cc] = nt.split()[1]
                chainDB[cc] = nt.split()[2]
            else:
                chainSeq[cc] += nt.split()[1]
                chainDB[cc] += nt.split()[2]
            missingchain = cc
            relabel[respos] = cc
        respos += 1
    chains = reorder
      
          
    # open the ct file and find (sub)chain interactions
    with open("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+id+".ct", 'r') as f:
        ctlines = f.readlines()
    for lp in range(1,len(ctlines)):
        lpbp = int(ctlines[lp].split()[4])
        if lpbp != 0:
            c1 = relabel[lp]
            c2 = relabel[lpbp]
            if c1 in chainBind:
                chainBind[c1].append(c2)
                if c1+'+'+c2 in BasePairs:
                    BasePairs[c1+'+'+c2].append(str(lp)+'+'+str(lpbp))
                else:
                    BasePairs[c1+'+'+c2] = [str(lp)+'+'+str(lpbp)]
            else:
                chainBind[c1] = [c2]
                BasePairs[c1+'+'+c2] = [str(lp)+'+'+str(lpbp)]
            if c2 in chainBind:
                chainBind[c2].append(c1)
                if c2+'+'+c1 in BasePairs:
                    BasePairs[c2+'+'+c1].append(str(lpbp)+'+'+str(lp))
                else:
                    BasePairs[c2+'+'+c1] = [str(lpbp)+'+'+str(lp)]
            else:
                chainBind[c2] = [c1]
                BasePairs[c2+'+'+c1] = [str(lpbp)+'+'+str(lp)]        


    # Reorganize chainBind dictionary
    for c in chainBind:        
        xRemove = []
      
        # Not include dimers
        cb = list(set(chainBind[c]))
        for x in cb:
            if x != c:
                if SequenceMatcher(None, chainSeq[x], chainSeq[c]).ratio() > 0.92:
                    if SequenceMatcher(None, chainDB[x], chainDB[c]).ratio() > 0.92:
                        xRemove.append(x)
                        Dimers.append(id)
      
        # Consider two chains interacting if bp>1
        cb = list(set(chainBind[c]))
        for x in cb:
            if c.split('-')[0] != x.split('-')[0]:
                if chainBind[c].count(x) <= 1:
                    xRemove.append(x)
      
        chainBind[c] = list(set(chainBind[c]))                        
                      
      
        # Need to erase these invalid or weak interactions
        xRemove = list(set(xRemove))
        for x in xRemove:
            chainBind[c].remove(x)
        for x in xRemove:
            BP = BasePairs[c+'+'+x]
            for bp in BP:
                p1 = int(bp.split('+')[0])
                p2 = int(bp.split('+')[1])
                # erase this base pair
                index = 0
                for i in range(len(db)):
                    if db[i] != '&':
                        index += 1
                        if index == p1 or index == p2:
                            db = list(db)
                            db[i] = '.'
                            db = ''.join(db)
                count = 0
                for ch in chains:
                    for i in range(len(chainDB[ch])):
                        count += 1
                        if count == p1 or count == p2:
                            chainDB[ch] = list(chainDB[ch])
                            chainDB[ch][i] = '.'    
                            chainDB[ch] = ''.join(chainDB[ch])


    PDBexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/' + id + '.pdb')
    if PDBexist:     
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb', 'r') as f:
            pdb = f.readlines()
    else:
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif', 'r') as f:
            pdb = f.readlines()
    i = 0
    t = []
    firsttitleline = True
    while i<len(pdb) and protein!='YES':
        l = pdb[i]
        i += 1
        if PDBexist:
            if l[0:5] == 'TITLE':
                if firsttitleline:
                    t += l.split()[1:]
                    firsttitleline = False
                else:
                    t += l.split()[2:]
        else:
            if 'title' in l and firsttitleline:
                t += l.split()[1:]
                firsttitleline = False
        if l[0:4] == 'ATOM':
            for aa in aaList:
                if aa in l:
                    protein = 'YES'
    for w in t:
        title = title + w + ' '

        
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'w') as f:
        f.write('Title\t\t\t'+title+'\n')
      
        s = ''
        for c in chains:
            s = s + c + ' '
        f.write('Chain order\t\t'+s+'\n')
      
        f.write('Amino acids\t\t'+protein+'\n')
      
        f.write('Metal')
        if metals:
            for me in metals:
                f.write('\t\t\t'+me+'\t'+metals[me]+'\n')
        else:
            f.write('\n')    
      
        f.write('Chains')
        for i in range(len(chains)):
            if chains[i] in chainType:
                f.write('\t\t\t'+chains[i]+'\t'+chainType[chains[i]]+'\t')
            else:
                f.write('\t\t\t'+chains[i]+'\t\t\t')
            ss = ''
            if chains[i] in chainBind:
                for cb in chainBind[chains[i]]:
                    ss = ss + cb + ' '
            f.write(ss+'\n')
      
        f.write('Pseudoknot\t\t'+pseudo+'\n')
      
        f.write('\nWhole structure\n')
        f.write(name)
        f.write(seq)
        f.write(db) 

        f.write('\nInteracting substructures\n')
      
        # Assign chain order number
        chainOrder = {}
        k = 0
        for c in chains:
            chainOrder[c] = k
            k += 1


      
################# New Algorithm to separate substructures, only include interacting chains ########################                
        lib = copy.copy(chains) # will remove a chain from library if assigned to a substructure
        subs = [] # record all the substructures
      
        for c in chains:
            if c in lib: # not assigned yet, find interacting chains
                ic = []
                if c in chainBind:
                    ic += chainBind[c]
                ic = list(set(ic))
                if c in ic:
                    ic.remove(c)
                # Find and include new chains that interact with the current substructure
                substructure = copy.copy(ic)
                new = copy.copy(ic)
                while new:
                    sub_new = copy.copy(substructure)
                    for x in new:
                        sub_new += chainBind[x]
                    sub_new = list(set(sub_new))
                    if c in sub_new:
                        sub_new.remove(c)
                    new = [s for s in sub_new if s not in substructure]
                    substructure = copy.copy(sub_new)
                substructure.append(c)
                substructure = sorted(substructure, key=lambda s: chainOrder[s])
                for x in substructure:
                    lib.remove(x)
                subs.append(substructure)
          
          
        # Write down the separated substructures
        for substructure in subs:
            ic = ''
            subseq = ''
            subdb = ''
            for sc in substructure:
                ic = ic + sc + ' '
                subseq += chainSeq[sc]
                subdb += chainDB[sc]
            f.write('Chain ' + ic + '\n')
            f.write(subseq + '\n')
            f.write(subdb + '\n\n')        
          
          
    # Write ct files for each substructure
    # os.system("export DATAPATH=/Users/qz886/Downloads/RNAstructure/data_tables")    
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
      
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            c = lines[i+1].split()[1:]
            seq = lines[i+2]
            db = lines[i+3]
            fname = id+'_'
            for x in c:
                fname += x
            with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/'+fname+'.out', 'w') as f:
                f.write('>'+fname+'\n')
                f.write(seq)
                f.write(db)
            os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct" )
            if not os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct"):
                Problems.append(fname)
            i += 5
          
            while i < len(lines):
                if lines[i] != '\n':
                    c = lines[i].split()[1:]
                    seq = lines[i+1]
                    db = lines[i+2]
                    fname = id+'_'
                    for x in c:
                        fname += x
                    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/'+fname+'.out', 'w') as f:
                        f.write('>'+fname+'\n')
                        f.write(seq)
                        f.write(db)
                    os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct" )
                    if not os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT/"+fname+".ct"):
                        Problems.append(fname)
                    i += 4
  
    print('Failed to produce ct files:')
    print(Problems)

print('Dimers:')
print(Dimers)



######################################################################################################################
########################## For large PDB files, use chains listed in BGSU file as filters ############################
######################################################################################################################
largePDB = []
largePDBfrag = []

with open('/Users/qz886/Desktop/Dual_Graph_Codes/Dual_Library_noGraph.txt', 'r') as f:
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



########################## Read the BGSU file to help organize large PDB files ##########################
with open('nrlist_3.209_all.csv', 'r') as f:
    lines = f.readlines()
   
BGSUchains = {}

for l in lines:
    name = l.split(',')[1]
    id = name.split('|')[0].split('\"')[1]
    c = []
    if '+' in name:
        allchains = name.split('+')
        for i in range(len(allchains)-1):
            c.append(allchains[i].split('|')[2])
        c.append(allchains[-1].split('|')[2].split('\"')[0])
    else:
        c.append(name.split('|')[2].split('\"')[0])
    if id in BGSUchains:
        BGSUchains[id].append(c)
    else:
        BGSUchains[id] = []
        BGSUchains[id].append(c)


         
for id in largePDB:
   
    print(id)

    name = ''
    seq = ''
    db = ''
    chains = []
    chainType = {}
    subchainPos = {} # for each chain, record which subchain every residue belongs to
    chainBind = {}
    BasePairs = {} # record base pairs between every two chains
    chainSeq = {}
    chainDB = {}
    metals = {}
    protein = 'NO'
    title = ''
    pseudo = 'NO'
    NT = [] # all nucleotide information
   
   
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR/'+id+'.out', 'r') as f:
        lines = f.readlines()
   
    n = 0    
    while n<len(lines):
       
        if 'no. of metals:' in lines[n]:
            me = lines[n].split()[-1]
            if me != '0':
                me = me.split('[')[1].split(']')[0].split(',')
                for m in me:
                    metals[m.split('=')[0]] = m.split('=')[1]
            n += 1

               
        elif lines[n] == 'Secondary structures in dot-bracket notation (dbn) as a whole and per chain\n':
            name = lines[n+1]
            seq = lines[n+2]
            db = lines[n+3]
            if '[' in db:
                pseudo = 'YES'
            n += 4            
            ############### New change, now consider broken subchains as individual chains ############################
            while lines[n][0] == '>':
                c = lines[n].split('>')[1].split(' ')[0].split('-')[1]
                subchains = lines[n+1].split('\n')[0].split('&')
                subchainDB = lines[n+2].split('\n')[0].split('&')
                if len(subchains) == 1:
                    chains.append(c)
                    chainType[c] = lines[n].split(' ')[-1].split('\n')[0]
                    subchainPos[c] = [c for s in subchains[0]]
                    chainSeq[c] = subchains[0]
                    chainDB[c] = subchainDB[0]
                else:
                    subchainPos[c] = []
                    for q in range(len(subchains)):
                        chains.append(c+'-'+str(q+1))
                        chainType[c+'-'+str(q+1)] = lines[n].split(' ')[-1].split('\n')[0]
                        subchainPos[c] += [c+'-'+str(q+1) for s in subchains[q]]
                        chainSeq[c+'-'+str(q+1)] = subchains[q]
                        chainDB[c+'-'+str(q+1)] = subchainDB[q]
                n += 3
               
        elif lines[n][0:30] == 'Summary of structural features':
            n += 8
            while lines[n] != '\n':
                NT.append(lines[n])
                n += 1
               
        else:
            n +=1


    # Write the original whole structure in ct format to help later chain reorganization
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/'+id+'.out', 'w') as f:
        f.write('>'+id+'\n')
        f.write(seq.replace('&',''))
        f.write(db.replace('&',''))
    os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+id+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+id+".ct" )
   

    #################### New change, now find all (sub)chain interactions using the ct file #############################
    # reorder the chains and relabel each residue because of subchains
    reorder = []
    relabel = {}
    missingchain = ''
    respos = 1
    for nt in NT:
        cc = nt.split()[3].split('.')[0]
        ################ New change, find DSSR files have errors, some chains are not listed in the Secondary structures section, e.g. 3AVT #############################
        if cc in subchainPos:
            subcc = subchainPos[cc].pop(0)
            if subcc not in reorder:
                reorder.append(subcc)
            relabel[respos] = subcc
        else:
            if cc != missingchain:
                reorder.append(cc)
                chainSeq[cc] = nt.split()[1]
                chainDB[cc] = nt.split()[2]
            else:
                chainSeq[cc] += nt.split()[1]
                chainDB[cc] += nt.split()[2]
            missingchain = cc
            relabel[respos] = cc
        respos += 1
    chains = reorder


          
    # open the ct file and find (sub)chain interactions
    with open("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+id+".ct", 'r') as f:
        ctlines = f.readlines()
    for lp in range(1,len(ctlines)):
        lpbp = int(ctlines[lp].split()[4])
        if lpbp != 0:
            c1 = relabel[lp]
            c2 = relabel[lpbp]
            if c1 in chainBind:
                chainBind[c1].append(c2)
                if c1+'+'+c2 in BasePairs:
                    BasePairs[c1+'+'+c2].append(str(lp)+'+'+str(lpbp))
                else:
                    BasePairs[c1+'+'+c2] = [str(lp)+'+'+str(lpbp)]
            else:
                chainBind[c1] = [c2]
                BasePairs[c1+'+'+c2] = [str(lp)+'+'+str(lpbp)]
            if c2 in chainBind:
                chainBind[c2].append(c1)
                if c2+'+'+c1 in BasePairs:
                    BasePairs[c2+'+'+c1].append(str(lpbp)+'+'+str(lp))
                else:
                    BasePairs[c2+'+'+c1] = [str(lpbp)+'+'+str(lp)]
            else:
                chainBind[c2] = [c1]
                BasePairs[c2+'+'+c1] = [str(lpbp)+'+'+str(lp)]        


    # Reorganize chainBind dictionary
    for c in chainBind:        
        xRemove = []
       
        # Not include dimers
        cb = list(set(chainBind[c]))
        for x in cb:
            if x != c:
                if SequenceMatcher(None, chainSeq[x], chainSeq[c]).ratio() > 0.92:
                    if SequenceMatcher(None, chainDB[x], chainDB[c]).ratio() > 0.92:
                        xRemove.append(x)
                        Dimers.append(id)
       
        # Consider two chains interacting if bp>1
        for x in cb:
            if c.split('-')[0] != x.split('-')[0]:
                if chainBind[c].count(x) <= 1:
                    xRemove.append(x)
       
        chainBind[c] = list(set(chainBind[c]))                        
                       
       
        # Need to erase these invalid or weak interactions
        xRemove = list(set(xRemove))
        for x in xRemove:
            chainBind[c].remove(x)
        for x in xRemove:
            BP = BasePairs[c+'+'+x]
            for bp in BP:
                p1 = int(bp.split('+')[0])
                p2 = int(bp.split('+')[1])
                # erase this base pair
                index = 0
                for i in range(len(db)):
                    if db[i] != '&':
                        index += 1
                        if index == p1 or index == p2:
                            db = list(db)
                            db[i] = '.'
                            db = ''.join(db)
                count = 0
                for ch in chains:
                    for i in range(len(chainDB[ch])):
                        count += 1
                        if count == p1 or count == p2:
                            chainDB[ch] = list(chainDB[ch])
                            chainDB[ch][i] = '.'    
                            chainDB[ch] = ''.join(chainDB[ch])       


    # Assign chain order number
    chainOrder = {}
    k = 0
    for c in chains:
        chainOrder[c] = k
        k += 1
         
    lib = copy.copy(chains) # will remove a chain from library if assigned to a substructure
    subs = [] # record all the substructures
   
    for c in chains:
        if c in lib: # not assigned yet, find interacting chains
            ic = []
            if c in chainBind:
                ic += chainBind[c]
            ic = list(set(ic))
            if c in ic:
                ic.remove(c)
            # Find and include new chains that interact with the current substructure
            substructure = copy.copy(ic)
            new = copy.copy(ic)
            while new:
                sub_new = copy.copy(substructure)
                for x in new:
                    sub_new += chainBind[x]
                sub_new = list(set(sub_new))
                if c in sub_new:
                    sub_new.remove(c)
                new = [s for s in sub_new if s not in substructure]
                substructure = copy.copy(sub_new)
            substructure.append(c)
            substructure = sorted(substructure, key=lambda s: chainOrder[s])
            for x in substructure:
                lib.remove(x)
            subs.append(substructure)
       
       
    # Only redo substructures too large in largePDB, use BGSU chains as filters
    redo_subs = []
   
    for substructure in subs:
        fname = id+'_'
        for x in substructure:
            fname += x

        # use BGSU filters to partition the subchains
        if fname in largePDBfrag:
            filters = []
            partitions = {}
            for BGSUsub in BGSUchains[id]:
                exist = False
                for sc in substructure:
                    if sc.split('-')[0] in BGSUsub:
                        exist = True
                        if str(BGSUsub) in partitions:
                            partitions[str(BGSUsub)].append(sc)
                        else:
                            partitions[str(BGSUsub)] = [sc]
                if exist:
                    filters.append(BGSUsub)
            # group all the subchains in each partition
            for fil in filters:
                par = partitions[str(fil)]
                # clean interactions for each subchain in this partition
                for sc in par:
                    xremove = []
                    for x in chainBind[sc]:
                        if x not in par:
                            xremove.append(x)
                    for x in xremove:
                        chainBind[sc].remove(x)
                        BP = BasePairs[sc+'+'+x]
                        for bp in BP:
                            p1 = int(bp.split('+')[0])
                            p2 = int(bp.split('+')[1])
                            index = 0
                            for i in range(len(db)):
                                if db[i] != '&':
                                    index += 1
                                    if index == p1 or index == p2:
                                        db = list(db)
                                        db[i] = '.'
                                        db = ''.join(db)
                            count = 0
                            for ch in chains:
                                for i in range(len(chainDB[ch])):
                                    count += 1
                                    if count == p1 or count == p2:
                                        chainDB[ch] = list(chainDB[ch])
                                        chainDB[ch][i] = '.'    
                                        chainDB[ch] = ''.join(chainDB[ch])
                # group into interacting substructures
                parOrder = {} # assign orders
                k = 0
                for c in par:
                    parOrder[c] = k
                    k += 1
                parlib = copy.copy(par) # will remove a chain from library if assigned to a substructure             
                for sc in par:
                    if sc in parlib: # not assigned yet, find interacting chains
                        ic = []
                        if sc in chainBind:
                            ic += chainBind[sc]
                        ic = list(set(ic))
                        if sc in ic:
                            ic.remove(sc)
                        # Find and include new chains that interact with the current substructure
                        parsub = copy.copy(ic)
                        new = copy.copy(ic)
                        while new:
                            sub_new = copy.copy(parsub)
                            for x in new:
                                sub_new += chainBind[x]
                            sub_new = list(set(sub_new))
                            if sc in sub_new:
                                sub_new.remove(sc)
                            new = [s for s in sub_new if s not in parsub]
                            parsub = copy.copy(sub_new)
                        parsub.append(sc)
                        parsub = sorted(parsub, key=lambda s: parOrder[s])
                        for x in parsub:
                            parlib.remove(x)
                        redo_subs.append(parsub)



    PDBexist = os.path.isfile('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/' + id + '.pdb')
    if PDBexist:     
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.pdb', 'r') as f:
            pdb = f.readlines()
    else:
        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_files/'+id+'.cif', 'r') as f:
            pdb = f.readlines()
    i = 0
    t = []
    firsttitleline = True
    while i<len(pdb) and protein!='YES':
        l = pdb[i]
        i += 1
        if PDBexist:
            if l[0:5] == 'TITLE':
                if firsttitleline:
                    t += l.split()[1:]
                    firsttitleline = False
                else:
                    t += l.split()[2:]
        else:
            if 'title' in l and firsttitleline:
                t += l.split()[1:]
                firsttitleline = False
        if l[0:4] == 'ATOM':
            for aa in aaList:
                if aa in l:
                    protein = 'YES'
    for w in t:
        title = title + w + ' '

          
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'w') as f:
        f.write('Title\t\t\t'+title+'\n')
       
        s = ''
        for c in chains:
            s = s + c + ' '
        f.write('Chain order\t\t'+s+'\n')
       
        f.write('Amino acids\t\t'+protein+'\n')
       
        f.write('Metal')
        if metals:
            for me in metals:
                f.write('\t\t\t'+me+'\t'+metals[me]+'\n')
        else:
            f.write('\n')    
       
        f.write('Chains')
        for i in range(len(chains)):
            if chains[i] in chainType:
                f.write('\t\t\t'+chains[i]+'\t'+chainType[chains[i]]+'\n')
            else:
                f.write('\t\t\t'+chains[i]+'\n')
       
        f.write('Pseudoknot\t\t'+pseudo+'\n')
       
        f.write('\nWhole structure\n')
        f.write(name)
        f.write(seq)
        f.write(db) 

        f.write('\nInteracting substructures\n')
        for substructure in redo_subs:
            ic = ''
            subseq = ''
            subdb = ''
            for sc in substructure:
                ic = ic + sc + ' '
                subseq += chainSeq[sc]
                subdb += chainDB[sc]
            f.write('Chain ' + ic + '\n')
            f.write(subseq + '\n')
            f.write(subdb + '\n\n') 

          
           
    # Write ct files for each substructure
    # os.system("export DATAPATH=/Users/qz886/Downloads/RNAstructure/data_tables/")    
    with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_2D_largePDB/'+id+'.txt', 'r') as f:     
        lines = f.readlines()
       
    for i in range(len(lines)):
        if lines[i] == 'Interacting substructures\n':
            if i+1 < len(lines):
                c = lines[i+1].split()[1:]
                seq = lines[i+2]
                db = lines[i+3]
                fname = id+'_'
                for x in c:
                    fname += x
                with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/'+fname+'.out', 'w') as f:
                    f.write('>'+fname+'\n')
                    f.write(seq)
                    f.write(db)
                os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct" )
                if not os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct"):
                    Problems.append(fname)
                i += 5
               
                while i < len(lines):
                    if lines[i] != '\n':
                        c = lines[i].split()[1:]
                        seq = lines[i+1]
                        db = lines[i+2]
                        fname = id+'_'
                        for x in c:
                            fname += x
                        with open('/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/'+fname+'.out', 'w') as f:
                            f.write('>'+fname+'\n')
                            f.write(seq)
                            f.write(db)
                        os.system("dot2ct /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".out /Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct" )
                        if not os.path.isfile("/Users/qz886/Desktop/Dual_Library_Update/PDB_DSSR_CT_largePDB/"+fname+".ct"):
                            Problems.append(fname)
                        i += 4
   
    print('Failed to produce ct files:')
    print(Problems)

print('Dimers:')
print(Dimers)