function Zv = Zv_update(C,Zv,Zd,Y,U,V,lata1,...
    beta1,miu1,Tv,Av,Iv,W)

beta11 = beta1.*(beta1>=0);  %beta1的正部
beta10 = (-beta1).*(beta1<0);   %beta1的负部

FZ = 2*(C.*C.*(Y.*W))*Zd*V*U'+lata1*Av*Zv*U*U'+beta11+miu1*Tv;
FM = 2*(C.*C.*((Zv*U*V'*Zd').*W))*Zd*V*U'+lata1*Iv*Zv*U*U'+beta10+miu1*Zv;
FM(FM<eps) = 10^-10;
Zv = FZ./FM.*Zv;  %更新Zv


end