function [U,V] = UV_update(C,Zv,Zd,Y,U,V,lata1,lata2,ar,Av,Ad,Iv,Id,W)

max_error = 10^-5;
for i=1:100
    U1 = U_update(C,Zv,Zd,Y,U,V,lata1,ar,Av,Iv,W);
    V1 = V_update(C,Zv,Zd,Y,U1,V,lata2,ar,Ad,Id,W);
    if norm(U-U1,'fro')/norm(U, 'fro')<max_error &...
            norm(V-V1,'fro')/norm(V, 'fro')<max_error
        break;
    end
    U = U1;
    V = V1;
end

end
