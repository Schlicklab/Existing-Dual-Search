This repository consists of codes and data for manuscript "RNA-As-Graphs Motif Atlas --- Dual Graph Library Update and Frameshifting Element Applications" by Qiyao Zhu and Tamar Schlick

To find existing dual graph and subgraph motifs from available RNA structures in Protein Data Bank:
1. Download BGSU representative PDB structure list: http://rna.bgsu.edu/rna3dhub/nrlist/. Here, we downloaded nrlist_3.209_all.csv.
2. Write the list of representative PDB IDs mentioned in the BGSU file, using script PDB_ID.py.
3. Download pdb/cif files of the representative PDBs, using the shell script batch_download.sh provided by PDB: https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script. Clean the downloaded files using PDB_ID.py.
4. Extract 2D structures (X3DNA-DSSR: https://x3dna.org/) from the pdb/cif files, using PDBto2D.py. Find independent substructures using the search algorithm in PDBto2D.py. This step gives results shared in folder PDB_DSSR_2D.
5. Identify dual graph motifs for the substructures, using Dual_Library.py.
6. For large substructures that have no dual graphs assigned (>9 helices), use BGSU Integrated Functional Elements (IFEs) as filters to do a second round of search. Find filtered substructures using PDBto2D.py. This step gives results shared in folders PDB_DSSR_2D_largePDB.
7. Identify dual graph motifs for the filtered substructures, using Dual_Library.py. Collect all information using Dual_Library_list.py. The results are written in Dual_Library.txt, which records all existing dual graph motifs with their weights.
8. Partition all substructures using Dual_Library_sub.py, and collect all existing dual subgraphs using Dual_Library_sub_list.py. The results are written in Dual_Library_Sub.txt.
