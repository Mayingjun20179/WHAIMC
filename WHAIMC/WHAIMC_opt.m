function P = WHAIMC_opt(W,train_mat, vir_sim,drug_sim,option)
Y = train_mat;
% c = option(1);  %
c = option(1);
C = ones(size(train_mat));
C(find(train_mat==1))=c;

num_factors = 50;  %
ar = option(2);  %
lata1 = option(3); %
lata2 = option(4); %

miu1 = 1;
miu2 = 1;


num_V = size(vir_sim,1);
num_D = size(drug_sim,1);
% 
num_K = 10;
[~,Av] = construct_Hypergraphs_knn(vir_sim,num_K);
[~,Ad] = construct_Hypergraphs_knn(drug_sim,num_K);
Av = (Av+Av')/2;
Ad = (Ad+Ad')/2;


% normalization
add_eps = 0.01;
Zv = diag(ones(num_V,1)./ sum(Av + add_eps,2))*(Av + add_eps);
Zd = (Ad + add_eps) *  diag(ones(1,num_D)./ sum(Ad + add_eps,1));

YY = Assist_method_Utils.WKNKN(Y, vir_sim,drug_sim,5, 0.8);
[U,V] = NNDSVD(YY,num_factors,0);
V = V';
% U = rand(num_V,num_factors);
% V = rand(num_D,num_factors);
Tv = Zv;
Td = Zd;
beta1 = zeros(size(Zv));
beta2 = zeros(size(Zd));
Iv = eye(size(Zv));
Id = eye(size(Zd));

%%%
gpu_transform;


numiter = 100;
%
loss = loss_function(Y,C,U,V,Zv,Zd,Iv,Id,Av,Ad,option,W);
loss_best = loss;
while numiter > 0
   
    %update U,V
    [U,V] = UV_update(C,Zv,Zd,Y,U,V,lata1,lata2,ar,Av,Ad,Iv,Id,W);


    % updateZv,Zd
    [Zv,Zd] = Zvd_update(C,Zv,Zd,Y,U,V,lata1,lata2,...
        beta1,beta2,miu1,miu2,Av,Ad,Iv,Id,Tv,Td,W);
    % updateTv,Td
    Tv = T_update(Zv, miu1,beta1);
    Td = T_update(Zd, miu2,beta2);
    % update beta1 and beta2
    beta1 = beta1 + miu1*(Tv-Zv);
    beta2 = beta2 + miu2*(Td-Zd);  
    
    %计算损失
    loss1 = loss_function(Y,C,U,V,Zv,Zd,Iv,Id,Av,Ad,option,W);
    loss = [loss;loss1];
    if abs(loss(end)-loss(end-1))<10^-6
        break;
    end
  
    
    if loss_best>loss1
        U_best = U;
        V_best = V;
        Zv_best = Zv;
        Zd_best = Zd;
        flag = 0;
        loss_best = loss1;
    else
        flag = flag+1;
        miu1 = miu1*0.8;
        miu2 = miu2*0.8;
    end
    numiter = numiter-1;
%     test_varialbe
end
% plot(log(loss))
P = Zv_best*U_best*V_best'*Zd_best';
P = gather(P);
end





















