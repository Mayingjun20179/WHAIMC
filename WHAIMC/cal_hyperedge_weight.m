%%%%计算超边权重
%S表示相似性矩阵
%H超图指示矩阵，行表示顶点，列表示边
function We = cal_hyperedge_weight(S,H)
%%%计算距离矩阵
S = (S+S')/2;
D = 1./(exp(S));
D(1:(size(D,1)+1):end) = 0;

[numV,nume] = size(H);  
wol = zeros(1,nume);
for j = 1:nume
    ind = find(H(:,j)>eps);  %第j条边对应的索引
    Dj = D(ind,ind);
    N = length(ind);
    P = zeros(N+1,N+1);
    P(1,:) = [0,ones(1,N)];
    P(:,1) = [0;ones(N,1)];
    P(2:end,2:end) = Dj;   %伪亲和矩阵
    wol(j) = sqrt(abs(det(P)))/(2^(N/2)*factorial(N));
end
miu = mean(wol);
We = exp(-wol/miu);

end