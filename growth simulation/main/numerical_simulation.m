function [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max)
% variables
% biomass, 1 for B.S, 2 for Nostoc
N1=[]; N2=[];
% c,n,p for the element
Rc=[];Rn=[];Rp=[];
% monitor birth rate (for f_i)
G1=[];G2=[];

% input: initial values

global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2

step=0.01;
time=0+step:step:t_max; 
for t=time
    % as step is small, we don't consider the update order
    % update f
    f1=min([MM(rc, Kc1, n1), MM(rn, Kn1, n1), MM(rp, Kp1, n1)]);
    % f1=min([MM_ori(rc, Kc1), MM_ori(rn, Kn1), MM_ori(rp, Kp1)]); 
    f2=min([MM(rn, Kn2, n2), MM(rp, Kp2, n2)]); % C has no limit on Nostoc
 
    % update nutrient concentration, and define Grow
    Grow1=f1*r1*n1;
    Grow2=f2*r2*n2;
%     Grow1=f1*r1*n1*toxin(ga1, n1, n2);
%     Grow2=f2*r2*n2*toxin(ga2, n1, n2);
    rc=rc+(-Re1*n1+NPh*n2-Q1c*Grow1-Q2c*(Grow2-m2*n2))*step;
    rn=rn+(Nf1*n1+Nf2*n2-Q1n*Grow1-Q2n*(Grow2-m2*n2))*step;
    rp=rp+(Pf*n1-Q1p*Grow1-Q2p*(Grow2-m2*n2))*step;
    % if dead B.S are decomposed 
%     rc=rc+(-Re1*n1+NPh*n2-Q1c*(Grow1-m1*n1)-Q2c*(Grow2-m2*n2))*step;
%     rn=rn+(Nf1*n1+Nf2*n2-Q1n*(Grow1-m1*n1)-Q2n*(Grow2-m2*n2))*step;
%     rp=rp+(Pf*n1-Q1p*(Grow1-m1*n1)-Q2p*(Grow2-m2*n2))*step;
       
    % update cell density
    n1=n1+(Grow1-m1*n1)*step;
    n2=n2+(Grow2-m2*n2)*step;
    
    % add values for main vars
    N1=[N1, n1];
    N2=[N2, n2];
    Rc=[Rc, rc];
    Rn=[Rn, rn];
    Rp=[Rp, rp];
    % monitored vars
    G1=[G1 Grow1];
    G2=[G2 Grow2];

end 
