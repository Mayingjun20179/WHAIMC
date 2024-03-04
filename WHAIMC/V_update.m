function V = V_update(C,Zv,Zd,Y,U,V,lata2,ar,Ad,Id,W)
CT = C';  YT = Y';
FZ = 2*Zd'*(CT.*CT.*(YT.*W'))*Zv*U+lata2*Zd'*Ad*Zd*V;
FM = 2*Zd'*(CT.*CT.*((Zd*V*U'*Zv').*W'))*Zv*U+ar*V+lata2*Zd'*Id*Zd*V;

FM(FM<eps) = 10^-10;
V = FZ./FM.*V;  %更新V

end