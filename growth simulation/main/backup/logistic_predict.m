function [u,C]=logistic_predict(x0,y,year_pre)
x1 = diff(x0);
x1=[x0(1),x1];
n = length(x0);           
z = (x0(2:n)+x0(1:n-1)) / 2;        %求x1的均值生成序列

B = [-z' z'.^2];
Y = x1(2:end)';
u = B\Y;                            %估计参数a,b的值

%创建符号变量
syms x(t)
x = dsolve(diff(x)+u(1)*x==u(2)*x^2, x(0)==x0(1));  
x=subs(x,{'a','b','x0'},{u(1),u(2),x1(1)});
x0_hat= subs(x,'t',[0:n-1+year_pre]);   %预测至2030年
x0_hat= double(x0_hat);

year_1 = y;
year_2 = [y,y(length(y))+1:y(length(y))+year_pre];

plot(year_1,x0)
hold on
plot(year_2,x0_hat)
grid on

epsilon=x0-x0_hat(1:n);
C=var(epsilon)/var(x0);
