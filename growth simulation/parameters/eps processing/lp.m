function [alpha,p0]=lp(A,Bb,P)
p0=P(1);
y=P-p0-Bb;
p=polyfit(A,y,1);
alpha=p(1);

% tests
% scatter(A,y)
% hold on
% plot(A,polyval(p,A))