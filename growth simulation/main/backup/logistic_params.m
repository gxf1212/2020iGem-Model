clear;clf;
data_population=xlsread('sign_data',5,'a30:b68');
Y=data_population(:,1);
P=data_population(:,2);
n=length(P)-1;
for t=1:n
	Z(t)=(P(t+1)-P(t))/P(t+1);
end
[B,Bint,r,rint,stats]=regress(Z',[ones(n,1) P(1:n)]);
gamma=B(1,1);
beta=B(2,1);
b=log(1-gamma);
c=-beta/gamma;
a=exp((sum(log(1./P(1:n)-c))-n*(n+1)*b/2)/n);
%forecast
y=[Y' Y(length(Y))+1:1:Y(length(Y))+30];
Pf=1./(c+a*exp(b*(y-y(1))));
plot(y,Pf,'r-o',Y(1:length(P)),P,'g-^');
