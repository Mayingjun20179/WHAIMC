function V_D = process_data( )
load('data\drug_virus.mat')
Nd = size(drug_virus.drug_sim.target_sim,1);
drug_virus.drug_sim.target_sim(1:(1+Nd):end) = 1;
V_D.drug_sim = {drug_virus.drug_sim.pubchem_sim,drug_virus.drug_sim.target_sim};
V_D.vir_sim = drug_virus.viru_sim.value;
%%%计算交互矩阵
[~,indd] = cellfun(@(x)ismember(x,drug_virus.drug_sim.name(:,1)),drug_virus.inter(:,1)); 
[~,indv] = cellfun(@(x)ismember(x,drug_virus.viru_sim.name(:,1)),drug_virus.inter(:,2)); 
inter_matrix = full(sparse(indv,indd,ones(length(indv),1)));
V_D.interaction = inter_matrix;  

end