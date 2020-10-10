global a21 a12 rho
a21=0.83;rho=1;a12=0.144; 
t1=0;
t2=20;
x0=0.1; %芽孢
y0=0.01;  %蓝藻
[t,z]=ode45(@myfun1,[t1,t2],[x0,y0]); 

subplot(2,1,1);
plot(t,z(:,1),t,z(:,2)); 
title('x,y关于t的函数图象'); 
subplot(2,1,2); 
plot(z(:,1),z(:,2));
title('x,y的相图');