function [L,AA]= construct_Hypergraphs_knn(W,k_nn)
%计算超图拉普拉斯和超图连接矩阵
% k_nn = 5;
L_H = [];
n_size = size(W,1);
n_vertex = n_size;
n_edge = n_size;
H = zeros(n_edge,n_vertex);
% W = (1/sqrt(2*pi))*exp(W);
%build Association matrix of Hypergraphs
for i=1:n_vertex
    ll = W(i,:);
    [B,index_i] = sort(ll,'descend');
    k_ii = index_i(1:k_nn);
%     H(i,k_ii) = 1;  %H的行表示超边，列表示顶点
    H(i,k_ii) = W(i,k_ii);  %H的行表示超边，列表示顶点
end

We = cal_hyperedge_weight(W,H');
We = diag(We);
% %%%%%方法1：计算Dv，De,
% H = H';
% [n,m] = size(H);
% Dv = H*diag(We); 
% DV_21 = zeros(n,1);
% DV_21(Dv>eps) = Dv(Dv>eps).^(-1/2);
% DV_21 = diag(DV_21);
% De = sum(H);
% De_1(De>eps) = De(De>eps).^(-1);
% De_1 = diag(De_1);
% AA = DV_21*H*We*De_1*H'*DV_21;
% L = eye(n)-AA;

%%%%%方法2：计算De,
H = H';
[n,m] = size(H);
De = sum(H);
De_1(De>eps) = De(De>eps).^(-1);
De_1 = diag(De_1);
AA = H*We*De_1*H';
AA = AA-diag(diag(AA));
Dn = diag(sum(AA));
Dn12 = Dn^(-1/2);
Dn12(isnan(Dn12)) = 0;
L = Dn12*(Dn - AA)*Dn12;
AA = Dn12*AA*Dn12;

end