% 模拟生长的总程序
% 使用真实参数
% 根据速率的单位，时间统一为min^{-1}
% 单位种群的速率、Ks也化为无量纲
% N的单位是g/L，虽然按照培养液的标准，但希望能假设土壤的菌和资源也用浓度

clear;clc;clf;
% 参数，不用调的设成全局变量
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf Pf Q1c Q1n Q1p Q2c Q2n Q2p

% 增长率，decay rate
r1=log(2)/30;
r2=log(2)/20;
m1=0.001;m2=0.001; 

% 半饱和值
Kc1=25;Kn1=9.31;Kp1=1.30; % 暂无单位种群的，单位g/L
Kn2=9.31;Kp2=1.30; % 单位为L^{-1}

% 转化速率，暂时认为和时间无关
% 芽孢呼吸，maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1=0.67e-3*12/60; 
% 蓝藻（用了铜绿）净光合，4.41+-0.16μmol(mg·min)，1000/54为光强
NPh=4.41e-3*12.01*tanh((1000/54)*0.016/4.41); 
% 蓝藻固氮,1.72+-0.25 fmol N/cell/h，细胞质量是估计的
Nf=1.72e-15*14.01/60*(7722e6/0.42); 
 Nf=Nf*10;
% 芽孢溶磷，无机+有机
Pf=8.87e-4/24/60+0; 
 Pf=Pf*10;

% Q值，假设为常数
%芽孢
Q1c=0.34*8/0.3*1e-3*72;
% Q1c=0.303; % 待定，能用
Q1n=0.0734;
Q1p=0.00607;
% 蓝藻
% Q2c=0.2523; 
Q2c=0.488/0.65;
Q2n=0.0846;
Q2p=0.00248;


%% 真实模拟
% 变量
% 生物量：1表示芽孢，2表示蓝藻
% c,n,p表示各自元素

% 迭代变量，设定为初值
% 初始菌浓度，通常藻比菌多
n1=0.01;n2=0.2;
% 初始营养浓度，单位g/L
% 源数据为有机质，非碳
% 沙丘
% rc=7.2e-3/1.724;rn=2.6e-4;rp=8.88e-7;
% 植被
rc=2.35e-2/1.724;rn=2.8e-4;rp=3.676e-6;
% 自己添加资源
% rc=3;rn=1;rp=0.1;

rho=1.3e3; %土壤密度，单位g/L
rc=rc*rho;rn=rn*rho;rp=rp*rho;


% 调用模拟函数
days=10;
t_max=1440*days; % 1天是1440min
[N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);

%% 展示结果
% 种群
figure(1)
plot(time, N1, time, N2);
legend('B.S','Nostoc')
% 资源
figure(2)
% plot(time, Rc, time, Rn, time, Rp);
% legend('C','N','P')
subplot(1,3,1)
plot(time, Rc);
legend('C')
subplot(1,3,2)
plot(time, Rn);
legend('N')
subplot(1,3,3)
plot(time, Rp);
legend('P')
% 生长、死亡速率对比
figure(3)
subplot(2,2,1)
plot(time, G1, time, N1*m1);
title('B.S');
legend('growth','decay')
subplot(2,2,2)
plot(time, G2, time, N2*m2);
title('Nostoc');
legend('growth','decay')
subplot(2,2,3)
plot(time, G1./N1/r1, time, G2./N2/r2);
title('f');
legend('B.S','Nostoc')


%% 调参


