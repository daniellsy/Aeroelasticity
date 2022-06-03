function [L_h,L_alpha,M_h,M_alpha]=parameter(k)
% 利用贝塞尔函数进行参数计算
i=sqrt(-1);
J0=besselj(0,k);
J1=besselj(1,k);
Y0=bessely(0,k);
Y1=bessely(1,k);
F1=J1*(J1+Y0)+Y1*(Y1-J0);
F2=(J1+Y0)^2+(Y1-J0)^2;
G1=-(Y1*Y0+J1*J0);
G2=(J1+Y0)^2+(Y1-J0)^2;
Fk=F1/F2;
Gk=G1/G2;
Ck=Fk+i*Gk;
L_h=1-2*i*Ck/k;
L_alpha=0.5-i*(1+2*Ck)/k-2*Ck/k/k;
M_h=0.5;
M_alpha=3/8-i/k;
end
