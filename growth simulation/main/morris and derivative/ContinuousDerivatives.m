% continuous derivatives with fixed delta mu
% original AdjustedMain_Derivatives

clear;
clc;
clf;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0
global n1 n2 rc rn rp t_max

% growth & decay rate
r1=log(2)/30;
r2=0.035/60;
m1=0.196/60; % 14.812 min^{-1} 
m2=0.0035/60;
r2=r2*1.13;

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
f_Nf=0.1;
f_Ph=1-f_Nf; 

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

% variables
% biomass, 1 for B.S, 2 for Nostoc
% c,n,p for the element

% initial values
% cell density, Nostoc > B.S often
n1=0.1;
n2=0.2;
% nutrient concentration, unit: g/L

% option 1: natural condition
% sand dune
rc=7.2e-3;rn=2.6e-4;rp=8.88e-7;rho=1.426e3;
% plants
% rc=2.35e-2;rn=2.8e-4;rp=3.676e-6;rho=1.446e3;

rc=rc/1.724; % I forget why?
rn=rn*0.05; % avalible N (inorganic and low-molecular-weight organic)
% rp = rp / 0.06  % just use soluble P

% rho: soil mass / volume in natural state. unit: g/L
% from percentage to mass / volume concentration. unit: g/L
rc=rc*rho;rn=rn*rho;rp=rp*rho; 

% option 2: trial concentration
% rc=0.02;rn=0.001;rp=0.001;rho=1.33e3;

% call the function
days=5;
t_max=1440*days; % 1 day=1440 min

%% simulation setup
% varsmax=31.0049/1440/1000; % 溶磷最大速率
% varsmin=varsmax*0.1; % 溶磷最小速率
% varsmax=varsmax*4;
% assign random values
% miu_num=100; % 共检测μ的数量
% times=1; % 每个μ正态随机取样点的个数
% vars: miu_num*times, storing all random values
% miu: 1*miu_num, storing all mu
% using min and max
% [vars,miu]=Norm(varsmax,varsmin,miu_num,times,-1,-1); % 获得所有的μ和对应的溶磷速率矩阵
% assign mu and sigma
% [vars,miu]=Norm(-1,-1,1,1,varsmax/2, varsmax/10);

varsmax=m1*1.2; % 溶磷最大速率
varsmin=m1*0.8; % 溶磷最小速率
% assign random values
miu_num=200; % 共检测μ的数量
times=1; % 每个μ正态随机取样点的个数
% vars: miu_num*times, storing all random values
% miu: 1*miu_num, storing all mu
% using min and max
% [vars,miu]=Norm(varsmax,varsmin,miu_num,times,-1,-1); % 获得所有的μ和对应的溶磷速率矩阵
% assign mu and sigma
% [vars,miu]=Norm(-1,-1,1,1,varsmax/2, varsmax/10);

% a linear vars, to see the parameters
vars=linspace(varsmin,varsmax,miu_num);
miu=linspace(varsmin,varsmax,miu_num);
delta_miu=varsmax/100000;

% 为了加快画图速度，也方便debug，特引入以下数组（用内存换时间）
t=zeros(1,miu_num*times); % 储存自变量
biomass_1=zeros(1,miu_num*times); % 储存第一个图的因变量
Biomass_1=zeros(1,miu_num*times); % 储存第一个图增加delta_miu的因变量
delta_1=zeros(1,miu_num*times); % 储存第一个图因变量的增加量
derivative_1=zeros(1,miu_num*times); % 储存第一个图的导数
biomass_2=zeros(1,miu_num*times); % 储存第二个图的因变量
Biomass_2=zeros(1,miu_num*times); % 储存第二个图增加delta_miu的因变量
delta_2=zeros(1,miu_num*times); % 存第二个图因变量的增加量
derivative_2=zeros(1,miu_num*times); % 储存第二个图的导数

for i=1:miu_num    
    t(1,i)=miu(i); % 储存对应的μ值
    [derivative_1(1,i),derivative_2(1,i),Biomass_1(1,i),Biomass_2(1,i),biomass_1(1,i),biomass_2(1,i),delta_1(1,i),delta_2(1,i)]=DerivativeContinuous(delta_miu,vars(i)); % 计算导数
end

%% 画图
figure(1)
subplot(1,2,1)
plot(t,biomass_1);
subplot(1,2,2)
scatter(t,derivative_1,'*');  % ,'size','0.5'

figure(2)
subplot(1,2,1)
plot(t,biomass_2);
subplot(1,2,2)
scatter(t,derivative_2,'*');

%% 标明图例
for i=1:2
    figure(i);
    subplot(1,2,1)
    xlabel('miu/(mg\cdot L^{-1}\cdot min^{-1})');
    ylabel('Biomass concentration/(g\cdot L^{-1})');
    if i==1
        title('B.S. Biomass Concentration');
    end
    if i==2
        title('Nostoc Biomass Concentration');
    end
    subplot(1,2,2)
    xlabel('miu/(mg\cdot L^{-1}\cdot min^{-1})');
    ylabel('Derivative of Biomass ');
    if i==1
        title('Derivative of B.S. Biomass Concentration');
    end
    if i==2
        title('Derivative of Nostoc Biomass Concentration');
    end
end


