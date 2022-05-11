# WHAIMC
Hypergraph Regularization and Adaptive Inductive Matrix Completion for Predicting Potential Drugs Against SARS-CoV-2
data: need to be decompressed
drug_name.xlsx: drug name
drug_pubchem_sim.txt: chemical structure similarity of drugs
drug_target_sim.txt: drug target similarity
inter.txt: virus-drug interaction matrix
virus_name.xlsx: virus name
virus_seq_similarity.txt: virus sequence similarity
drug_virus.mat: Aggregated matlab data file


code:
WHAINMF.zip: WHAIMC related files, need to be decompressed
mianZ.m: main function for cross-experiment
Before running the program, you need to download "NMFLibrary-master", put it in the current running path, import all files into the matlab path, run mianZ, and you can get the corresponding results.
