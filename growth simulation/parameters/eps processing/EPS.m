clear;clf;clc;

% data
datab=[ 1.13378684807257e-002    5.83879761290798e-003
 5.00000000000000e+000    1.75384955426214e-002
 1.00028344671202e+001    3.50106829735504e-002
 1.50056689342404e+001    7.26994769026744e-002
 2.00085034013605e+001    7.35651661386576e-002
 2.49971655328798e+001    8.30987990864216e-002
 2.99716553287982e+001    1.02744419067266e-001];
datae=[ 1.13378684807257e-002    5.83879761290798e-003
 5.00000000000000e+000    2.25926471671701e-002
 1.00170068027211e+001    4.13346349370073e-001
 1.49348072562358e+001    6.75605982465188e-001
 1.99943310657596e+001    6.78623001547189e-001
 2.49688208616780e+001    9.69748765932366e-001
 2.99433106575964e+001    1.59589258085906e+000];
datab=datab(1:5,:); % select the first five points
datae=datae(1:5,:);

len=length(datab(:,1)); % num of columns
t=datab(:,1);
X=datab(:,2);
P=datae(:,2);

% fitting
[mu,xm,x0]=growth(X,t);
dpdt_stable=(datae(len,2)-datae(len-1,2))/(datae(len,1)-datae(len-1,1));
beta=dpdt_stable/xm;

A=X-x0;
Bb=beta*xm/mu*log(1-xm/x0*(1-exp(mu*t)));
[alpha,p0]=lp(A,Bb,P);

%% final test
t_fit=0:0.2:30;
eT=exp(mu*t_fit);
X_fit=x0*eT./(1-x0/xm*(1-eT));
P_fit=p0+alpha*(X_fit-x0)+beta*xm/mu*log(x0*eT./X_fit);

clf;figure(1)
plot(t_fit,X_fit,t_fit,P_fit,'linewidth',1.5)
hold on
scatter(t,X,'*')
scatter(t,P,'*')
title('eps data fitting');
legend('biomass','EPS','biomass','EPS')
xlabel('time/days');
ylabel('concentration/(g\cdot L^{-1})');
