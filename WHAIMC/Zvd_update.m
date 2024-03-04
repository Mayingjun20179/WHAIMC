function [Zv,Zd] = Zvd_update(C,Zv,Zd,Y,U,V,lata1,lata2,...
    beta1,beta2,miu1,miu2,Av,Ad,Iv,Id,Tv,Td,W)

max_error = 10^-3;
for i=1:100
    Zv1 = Zv_update(C,Zv,Zd,Y,U,V,lata1,beta1,miu1,Tv,Av,Iv,W);
    Zd1 = Zd_update(C,Zv1,Zd,Y,U,V,lata2,beta2,miu2,Td,Ad,Id,W);
    
    if norm(Zv-Zv1,'fro')/norm(Zv,'fro')<max_error &...
            norm(Zd-Zd1,'fro')/norm(Zd,'fro')<max_error
        break;
    end
    Zv = Zv1;
    Zd = Zd1;
    
end


end
