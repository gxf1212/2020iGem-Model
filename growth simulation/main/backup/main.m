global a21 a12 rho
a21=0.83;rho=1;a12=0.144; 
t1=0;
t2=20;
x0=0.1; %ѿ��
y0=0.01;  %����
[t,z]=ode45(@myfun1,[t1,t2],[x0,y0]); 

subplot(2,1,1);
plot(t,z(:,1),t,z(:,2)); 
title('x,y����t�ĺ���ͼ��'); 
subplot(2,1,2); 
plot(z(:,1),z(:,2));
title('x,y����ͼ');