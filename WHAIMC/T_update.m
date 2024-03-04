function [T] = T_update(Z, miu,beta)

Z = gather(Z);
beta = gather(beta);
%%%%¸üÐÂTvºÍTd
[n,m] = size(Z);

d = Z-beta/miu;

T = zeros(n,m);

for i = 1:n
    dv = d(i, :);
    %     S(i, :) = EProjSimplex_new(dv);
    T(i, :) = rojSimplex_H(dv); % my own
end

T = gpuArray(single(T));
end