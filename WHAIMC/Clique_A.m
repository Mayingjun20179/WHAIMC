function C = Clique_A(W,H)
%W 表示边的权重
%H 行表示边 列表示点
%N 表示边的个数
N = size(H,1);  %边的个数
A=eye(N);%%A为对角线是1的对角矩阵
if (size(W,2)==1)
    W=diag(W);
end
%X=1/3*H'*W*H;
X=H'*W*(inv(2*A))*H;%%X=HW1/2H因为是3均匀图，所以De是3减去单位矩阵变为2
C = X-diag(diag(X));
C = sparse(C);
end