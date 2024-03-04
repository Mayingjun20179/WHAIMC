%%%将本研究的数据转化为GCN需要的数据
clc
clear
load('drug_virus.mat')
drug_name = drug_virus.drug_sim.name(:,1);
virus_name = drug_virus.viru_sim.name(:,1);
[~,indd] = cellfun(@(x)ismember(x,drug_virus.drug_sim.name(:,1)),drug_virus.inter(:,1)); 
[~,indv] = cellfun(@(x)ismember(x,drug_virus.viru_sim.name(:,1)),drug_virus.inter(:,2)); 
inter = full(sparse(indv,indd,ones(length(indv),1)));


fid = fopen('.\inter.txt','wt');     % 要输出文本的路径位置及名称
[m,n] = size(inter);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%d',inter(i,j));
        if j<n
            fprintf(fid,',');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

drug_sim = drug_virus.drug_sim.pubchem_sim;
fid = fopen('.\drug_pubchem_sim.txt','wt');     % 要输出文本的路径位置及名称
[m,n] = size(drug_sim);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f',drug_sim(i,j));
        if j<n
            fprintf(fid,',');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

drug_sim = drug_virus.drug_sim.target_sim;
fid = fopen('.\drug_target_sim.txt','wt');     % 要输出文本的路径位置及名称
[m,n] = size(drug_sim);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f',drug_sim(i,j));
        if j<n
            fprintf(fid,',');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);



virus_sim = drug_virus.viru_sim.value;
fid = fopen('.\virus_seq_similarity.txt','wt');     % 要输出文本的路径位置及名称
[m,n] = size(virus_sim);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%f',virus_sim(i,j));
        if j<n
            fprintf(fid,',');
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

Nd = length(drug_name);
drug_name = [num2cell([1:Nd]'),drug_name];
drug_name = [{'ID','Drugs'};drug_name];
writecell(drug_name,'.\drug_name.xlsx');

Nv = length(virus_name);
virus_name = [num2cell([1:Nv]'),virus_name];
virus_name = [{'ID','Virus'};virus_name];
writecell(virus_name,'.\virus_name.xlsx');