function [mu,xm,x0]=growth(X,t)

% T=exp(t);
% x1=1./X;
% t1=1./T;
% p=polyfit(t1,x1,1);
% xm=1/p(2);
x0=X(1);
% mu=log((1/x0-1/xm)/p(1));

% xm=0.5264; 
xm=0.10376; % a tried self-consistant value
p=polyfit(t,log(X./(xm-X)),1);
mu=p(1);

% tests
% log(xm/x0-1)+p(2)
% scatter(t,log(X./(xm-X)))





