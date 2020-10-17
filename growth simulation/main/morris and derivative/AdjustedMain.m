% main for growth simulation
% real parameters
% unit of time is min^{-1} in accord with those of rates
% rates per biomass, Ks are deminsionless（无量纲）
% we assume that the unit of N is g/L

% MATLAB supports versions Python 2.7, 3.6, and 3.7

clear;clc;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0

% growth & decay rate
r1=log(2)/30;
r2=0.03/60;
m1=0.196/60; % 14.812 min^{-1} 
% m1=0.003; 
m2=0.0035/60;
% m2=0.00;
% r2=r2*0.34;
% m2=m2*0.34;

% toxinic supression, gamma, unit: (g/L)^{-2}
ga1=0.1; % 2 on 1
ga2=0.1;

% half-saturation constants, unit: L^{-1}
% B. S
% Kc1=7.4482*0.4; % batch
% Kc1=4.3388*0.4; % fed-batch
Kc1=2.865*0.4;
Kn1=0.638;
Kp1=0.057;
% cyanobacteria
% Kn2=9.31;
Kn2=0.050993; % 蓝丝菌lansijun
Kp2=0.018739; % 铜绿的tonglvjiadanbao

% Q, constant
% B. S
% Q1c=0.303; % estimated elemental mass fraction
Q1c=0.34*8/0.3*1e-3*72; % reciprocal of growth/biomass yield
Q1n=0.0734;
Q1p=0.00607;
% Nostoc
% Q2c=0.2523; % estimated elemental mass fraction
Q2c=0.488/0.65; % not sensitivity, little change in result
Q2n=0.0846;
Q2p=0.00248;
% Q2p=Q2p*20; % consider P luxury uptake by Nostoc 
% Q1c=Q1c*20;

% To be realistic, we divide Nostoc into vegetative and N-fixing cells
f_Nf=0.1;f_Ph=1-f_Nf; 

% rates, as constant
% B. S respiration, maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1=0.67e-3*12/60; 
% (pseudomonas aeruginosa) net photosynthesis rate, 4.41+-0.16μmol(mg·min)
I=100000/54; % luminosity
% NPh=4.41e-3*12.01*tanh(I*0.016/4.41); % 0.0530
NPh=4.41e-3*12.01*I/(165+I+I^2/457); % 0.0103
NPh=NPh*f_Ph;
% bacteria, N decomposition
Nf1=2.600e-4;
% Nf1=2.3e-4;
% Nf1=2.202e-4;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/h，estimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
% B.S ^{-1}, inorganic+organic
% Pf=8.87e-4/1440+0; 
% Pf=Pf*10; 
Pfmax=31.0049/1440/1000; % 溶磷最大速率
Pfmin=Pfmax*0.1; % 溶磷最小速率
miu_num=100; % 共检测μ的数量
times=10; % 每个μ正态随机取样点的个数
[PF,miu]=Norm(Pfmax,Pfmin,miu_num,times,-1,-1); % 获得所有的μ和对应的溶磷速率矩阵

% pH
pH0=9.1;

%% simulation setup
% variables
% biomass, 1 for B.S, 2 for Nostoc
% c,n,p for the element

% initial values
% cell density, Nostoc > B.S often
n1=0.01;
n2=0.02;
% nutrient concentration, unit: g/L

% option 1: natural condition
% sand dune
rc=7.2e-3;rn=2.6e-4;rp=8.88e-7;rho=1.426e3;
% plants
% rc=2.35e-2;rn=2.8e-4;rp=3.676e-6;rho=1.446e3;

% rc=rc/1.724; % I forget why?
rn=rn*0.05; % avalible N (inorganic and low-molecular-weight organic)
% rp = rp / 0.06  % just use soluble P

% rho: soil mass / volume in natural state. unit: g/L
% from percentage to mass / volume concentration. unit: g/L
rc=rc*rho;rn=rn*rho;rp=rp*rho; 

% option 2: trial concentration
% rc=0.02;rn=0.001;rp=0.001;rho=1.33e3;

% option 3: cultural medium * 2
% 1 and 2 now means type of medium culture
% % concentration or percentage
% cm1=1;cm2=1;
% % mass fraction
% rc1=1;rn1=1;rp1=1;
% rc2=1;rn2=1;rp2=1;
% % mixture
% rc=cm1*rc1+cm2+rc2;
% rn=cm1*rn1+cm2+rn2;
% rp=cm1*rp1+cm2+rp2;

% call the function
days=5;
t_max=1440*days; % 1 day=1440 min
% 为了加快画图速度，特引入以下三个数组
t=zeros(1,miu_num*times); % 储存自变量
biomass_1=zeros(1,miu_num*times); % 储存第一个图的因变量
biomass_2=zeros(1,miu_num*times); % 储存第二个图的因变量
for i=1:miu_num
    color1=(1/miu_num)*i; % 设定渐变的颜色
    color2=(1/miu_num)*(miu_num-i);
    color3=0;
    color=[color1,color2,color3];
    for j=1:times
        Pf=PF(i,j); % 在每一个溶磷速率下进行一次计算
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        x=ones(1,length(time()))*miu(i); % 用μ的值作为x轴
%         figure(1); % 作芽孢杆菌的图像
%         plot3(x,time, N1, 'linewidth', 0.5,'color',color);
%         hold on
%         figure(2); % 作蓝藻的图像
%         plot3(x,time, N2, 'linewidth', 0.5,'color',color);
%         hold on
        biomass_1(1,(i-1)*j+j)=N1(length(N1)); % 储存第五天的生物量
        biomass_2(1,(i-1)*j+j)=N2(length(N2));
        t(1,(i-1)*j+j)=miu(i); % 储存对应的μ值
    end
end

% 画图
figure(1)
semilogy(t,biomass_1,'*','markersize',3);
hold on;
figure(2)
semilogy(t,biomass_2,'*','markersize',3);
hold on;

% 标明图例
for i=1:2
      figure(i);
%     title('growth curve');
%     ylabel('time/min');
%     xlabel('miu/mg*L^{-1}*min^{-1}');
%     zlabel('biomass concentration/(g\cdot L^{-1})');
%     if i==1
%         title('B.S. Growth Curve');
%     else
%         title('Nostoc Growth Curve');
%     end
    xlabel('miu/(mg\cdot L^{-1}\cdot min^{-1})');
    ylabel('biomass concentration/(g\cdot L^{-1})');
    if i==1
        title('B.S. Biomass Concentration');
    else
        title('Nostoc Biomass Concentration');
    end
end


