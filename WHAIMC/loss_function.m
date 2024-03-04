function loss = loss_function(Y,C,U,V,Zv,Zd,Iv,Id,Av,Ad,option,W)

ar = option(2);  %特征正则化参数
lata1 = option(3); %病毒正则化参数
lata2 = option(4); %药物正则化参数
%%%%%目标1
f1 = C.*(Y-Zv*U*V'*Zd').*W;
f1 = norm(f1,'fro')^2;

%%%%%目标2
f2 = ar/2 * (norm(U,'fro')^2+norm(V,'fro')^2);

%%%%%目标3
Lv = Iv - Av;
Ld = Id - Ad;
f3 = lata1 * trace((Zv*U)'*Lv*Zv*U)+...
    lata2 * trace((Zd'*V)'*Ld*Zd*V);

%%%%%总目标
loss = f1+f2+f3;

end