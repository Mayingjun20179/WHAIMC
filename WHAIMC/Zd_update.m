function Zd = Zd_update(C,Zv,Zd,Y,U,V,lata2,...
    beta2,miu2,Td,Ad,Id,W)
beta21 = beta2.*(beta2>=0);  %beta2的正部
beta20 = (-beta2).*(beta2<0);   %beta2的负部

FZ = 2*(C'.*C'.*(Y'.*W'))*Zv*U*V'+lata2*Ad*Zd*V*V'+beta21+miu2*Td;
FM = 2*(C'.*C'.*((Zd*V*U'*Zv').*W'))*Zv*U*V'+lata2*Id*Zd*V*V'+beta20+miu2*Zd;
FM(FM<eps) = 10^-10;
Zd = FZ./FM.*Zd;  %更新Zv


end