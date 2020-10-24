% main for growth simulation
% real parameters
% unit of time is min^{-1} in accord with those of rates
% rates per biomass, Ks are deminsionless（无量纲）
% we assume that the unit of N is g/L

% MATLAB supports versions Python 2.7, 3.6, and 3.7

clear;
clc;
clf;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0
global n1 n2 rc rn rp t_max step

% growth & decay rate
r1=log(2)/30;
r2=0.035/60;
m1=0.196/60; % 14.812 min^{-1} 
m2=0.0035/60;
r2=r2*1.13;

% toxinic supression, gamma, unit: (g/L)^{-2}
ga1=0.1; % 2 on 1
ga2=0.1;

% half-saturation constants, unit: L^{-1}
% B. S
Kc1=2.865*0.4;
Kn1=0.638;
Kp1=0.057;
% cyanobacteria
Kn2=0.050993; % lansijun
Kp2=0.018739; % tonglvjiadanbao

% Q, constant
% B. S
Q1c=0.34*8/0.3*1e-3*72; % reciprocal of growth/biomass yield
Q1n=0.0734;
Q1p=0.00607;
% Nostoc
Q2c=0.488/0.65; % not sensitivity, little change in result
Q2n=0.0846;
Q2p=0.00248;

% To be realistic, we divide Nostoc into vegetative and N-fixing cells
f_Nf=0.1;
f_Ph=1-f_Nf; 

% rates, as constant
% B. S respiration, maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1=0.67e-3*12/60; 
% (pseudomonas aeruginosa) net photosynthesis rate, 4.41+-0.16μmol(mg·min)
I=100000/54; % luminosity
NPh=4.41e-3*12.01*I/(165+I+I^2/457); % 0.0103
NPh=NPh*f_Ph;
% bacteria, N decomposition
Nf1=2.600e-4;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/h，estimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
% B.S ^{-1}, inorganic+organic

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
rc=7.2e-3;
rn=2.6e-4;
rp=8.88e-7;
rho=1.426e3;
% plants

rn=rn*0.05; % avalible N (inorganic and low-molecular-weight organic)

% rho: soil mass / volume in natural state. unit: g/L
% from percentage to mass / volume concentration. unit: g/L
rc=rc*rho;
rn=rn*rho;
rp=rp*rho;

% call the function
days=6;
t_max=1440*days; % 1 day=1440 min

% simulation setup
Pf=31.0049/1440/1000; % 溶磷最大速率
% assign random values
miu_num=1000; % 共检测μ的数量
times=1; % 每个μ正态随机取样点的个数
step=50;
as1="B.S.";
as2="Nostoc";

% a linear Pf, to see the parameters

name="Nf1";

if strcmp(name,"Nf1")==1
    Max=Nf1*1.5;
    min=Nf1*0.01;
end

if strcmp(name,"Pf")==1
    Max=Pf*4;
    min=Pf*0.1;
end

VAR=linspace(min,Max,miu_num);
miu=linspace(min,Max,miu_num);
delta_miu=(Max-min)/1000000000000;

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
    [derivative_1(1,i),derivative_2(1,i),Biomass_1(1,i),Biomass_2(1,i),biomass_1(1,i),biomass_2(1,i),delta_1(1,i),delta_2(1,i)]=Derivative(delta_miu,VAR(i),name); % 计算导数
end

% 画图&标图例&调整

if strcmp(name,"Pf")==1 % Changes of Biomass Affected by Phosphorus Dissolving Rate
    
    DrawFigures("B.S. Biomass Concentration",t,biomass_1,"Phosphorus Dissolving Rate/(mg\cdot L^{-1}\cdot min^{-1})","Biomass Concentration/(g\cdot L^{-1})",as1,0);
    export_fig B.S._Biomass_Concentration -transparent;
    export_fig B.S._Biomass_Concentration.svg -transparent;
    
    DrawFigures("Derivative of B.S. Biomass Concentration with Phosphorus Dissolving Rate",t,derivative_1,"Phosphorus Dissolving Rate/(mg\cdot L^{-1}\cdot min^{-1})","Derivative of Biomass Concentration/(min^{-1})",as1,-inf);
    export_fig Derivative_of_B.S._Biomass Concentration_with_Phosphorus_Dissolving_Rate -transparent;
    export_fig Derivative_of_B.S._Biomass Concentration_with_Phosphorus_Dissolving_Rate.svg -transparent;
    
    DrawFigures("Nostoc Biomass Concentration",t,biomass_2,"Phosphorus Dissolving Rate/(mg\cdot L^{-1}\cdot min^{-1})","Biomass Concentration/(g\cdot L^{-1})",as2,0);
    export_fig Nostoc_Biomass_Concentration -transparent;
    export_fig Nostoc_Biomass_Concentration.svg -transparent;
    
    DrawFigures("Derivative of Nostoc Biomass Concentration with Phosphorus Dissolving Rate",t,derivative_2,"Phosphorus Dissolving Rate/(mg\cdot L^{-1}\cdot min^{-1})","Derivative of Biomass Concentration/(min^{-1})",as2,-inf);
    export_fig Derivative_of_Nostoc_Biomass_Concentration_with_Phosphorus_Dissolving_Rate -transparent;
    export_fig Derivative_of_Nostoc_Biomass_Concentration_with_Phosphorus_Dissolving_Rate.svg -transparent;
    
end

if strcmp(name,"Nf1")==1 % Changes of Biomass Affected by B.S. Nitrogen Decomposition Rate
    
    DrawFigures("B.S. Biomass Concentration",t,biomass_1,"Nitrogen Decomposition Rate/(mg\cdot L^{-1}\cdot min^{-1})","Biomass Concentration/(g\cdot L^{-1})",as1,0);
    export_fig B.S._Biomass_Concentration -transparent;
    export_fig B.S._Biomass_Concentration.svg -transparent;
    
    DrawFigures("Derivative of B.S. Biomass Concentration with B.S. Nitrogen Decomposition Rate",t,derivative_1,"Nitrogen Decomposition Rate/(mg\cdot L^{-1}\cdot min^{-1})","Derivative of Biomass Concentration/(min^{-1})",as1,-inf);
    export_fig Derivative_of_B.S._Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate -transparent;
    export_fig Derivative_of_B.S._Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate.svg -transparent;
    
    DrawFigures("Nostoc Biomass Concentration",t,biomass_2,"Nitrogen Decomposition Rate/(mg\cdot L^{-1}\cdot min^{-1})","Biomass Concentration/(g\cdot L^{-1})",as2,0);
    export_fig Nostoc_Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate -transparent;
    export_fig Nostoc_Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate.svg -transparent;
    
    DrawFigures("Derivative of Nostoc Biomass Concentration with B.S. Nitrogen Decomposition Rate",t,derivative_2,"Nitrogen Decomposition Rate/(mg\cdot L^{-1}\cdot min^{-1})","Derivative of Biomass Concentration/(min^{-1})",as2,-inf);
    export_fig Derivative_of_Nostoc_Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate -transparent;
    export_fig Derivative_of_Nostoc_Biomass_Concentration_with_B.S._Nitrogen_Decomposition_Rate.svg -transparent;
        
end
    


