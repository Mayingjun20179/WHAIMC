function U = U_update(C,Zv,Zd,Y,U,V,lata1,ar,Av,Iv,W)

FZ = 2*Zv'*(C.*C.*(Y.*W))*Zd*V+lata1*Zv'*Av*Zv*U;
FM = 2*Zv'*(C.*C.*((Zv*U*V'*Zd').*W))*Zd*V+ar*U+lata1*Zv'*Iv*Zv*U;

FM(FM<eps) = 10^-10;
U = FZ./FM.*U;  %更新U

end
