This Matlab code is from "Yingjun Ma, Junjiang Zhong, Nenghui Zhu. Weighted hypergraph learning and adaptive inductive matrix completion for SARS-CoV-2 drug repositioning  [J]. Methods, 2023, (219), 102-110. The Doi: 10.1016 / j.y meth. 2023.10.002. ", 
if the reference the code, please refer to the corresponding literature.
(Written by Yingjun Ma 2023)


To run the code:
1. Change Matlab work directory to "/WHAIMC/".

2. Run  "loadpah" code to add the current folder and subfolders into Matlab path searching list.

3. Open and run the demo file. 


Demo code:   mianZ.m

The "main_cv.m" contains three 5-CV scenarios: "Pairwise association" scenario,
 "New viruses" scenario, and "New drugs" scenario.


The “WHAIMC” folder contains the relevant calculation code for WHAIMC:

WHAIMC_opt：Weighted hypergraph learning and adaptive inductive matrix completion.

The Dataset folder contains experimental data in this paper:

(1)drug_name.xlsx：the names of all the drugs

(2)drug_pubchem_sim.txt: chemical molecular similarity of drugs

(3)drug_target_sim.txt:  drug similarity calculated according to drug-target association and target similarity

(4)inter.txt: drug-virus associations

(5)virus_name.xlsx: the names of all the virus

(6) virus_seq_similarity.txt: the virus similarity was calculated based on the complete genome sequences of the virus

(7)drug_virus.mat:summary matlab data file, including drug-virus association, drug similarity, and virus similarity




In this package, we used the NMFLibrary, which is downloaded from (https://github.com/hiroyuki-kasai/NMFLibrary/tree/master)
