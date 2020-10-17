% calculate derivative with multiple mu
% to find relationships between parameters?

clear;clc;clf;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0

% growth & decay rate
r1=log(2)/30;
r2=0.035/60;
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
Pf=31.0049/1440/1000; 

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

% rc=rc/1.724; % I forget why?
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

%% simulation of Nf1
% % variable setup
% Nf1=2.600e-4;
% 
% % task 1: fix mu, find the distribution of derivatives
% % miu: miu_num*1, storing all mu
% miu_num=7; % evenly spaced
% miu=linspace(Nf1*0.95,Nf1*1.05,miu_num)';
% % miu=Nf1; % manual
% % miu_num=length(miu); 
% 
% sigma=Nf1*0.01;
% % step=Nf1*0.005;
% % 每个μ正态随机取样点的个数
% times=40; 
% 
% vars=Norm_deri(miu,times,sigma,-1); % mode1, single sigma
% % vars=Norm_deri(miu,times,-1,); % step
% 
% % each row for a given miu and its different "times", or delta mu
% derivative_1=zeros(miu_num,times); 
% derivative_2=zeros(miu_num,times); 
% 
% for i=1:miu_num
%     Nf1=miu(i);
%     [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
%     t_total=length(time);
%     N1_cen=N1(t_total);
%     N2_cen=N2(t_total);
%     for j=1:times
%         Nf1=vars(i,j); % 在每一个溶磷速率下进行一次计算
%         [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
%         derivative_1(i,j)=(N1(t_total)-N1_cen)/(vars(i,j)-miu(i)); 
%         derivative_2(i,j)=(N2(t_total)-N2_cen)/(vars(i,j)-miu(i));
%     end
% end

% pause;

%% simulation of Pf
Pfmax=31.0049/1440/1000; % 溶磷最大速率
Pfmin=Pfmax*0.1; % 溶磷最小速率

% task 1: fix mu, find the distribution of derivatives
% miu: miu_num*1, storing all mu
miu_num=5; 
miu=linspace(Pfmax/4,Pfmax/4*5,miu_num)';
% 每个μ正态随机取样点的个数
times=40; 
sigma=Pf*0.01; % you should pick a low sigma if you want to see derivatives

% vars: miu_num*times, storing all random values (delta mu)
vars=Norm_deri(miu,times,Pfmax*0.2,-1); % mode1, single sigma
% vars=Norm_deri(miu,times,-1,Pfmax*0.01); % step

% each row for a given miu and its different "times", or delta mu
derivative_1=zeros(miu_num,times); 
derivative_2=zeros(miu_num,times); 

for i=1:miu_num
    Pf=miu(i);
    [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
    t_total=length(time);
    N1_cen=N1(t_total);
    N2_cen=N2(t_total);
    for j=1:times
        Pf=vars(i,j); % 在每一个溶磷速率下进行一次计算
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        derivative_1(i,j)=(N1(t_total)-N1_cen)/(Pf-miu(i)); 
        derivative_2(i,j)=(N2(t_total)-N2_cen)/(Pf-miu(i));
    end
end

% pause;

%% 画图分析
% mode 1
% for each mu, compare their derivatives
% figure(1)
clc;
for i=1:miu_num
    data1=derivative_1(i,:);
    data2=derivative_2(i,:);
    % plot
    figure(i)
    clf;
    hold on
    scatter(vars(i,:)-miu(i),data1,'*');
    scatter(vars(i,:)-miu(i),data2,'o');
    legend("B.S","Nostoc");
    xlabel('delta var');
    ylabel('derivative')
    miui=string(miu(i)*10^4)+'* 10^{-4}';
    title("var="+miui); % Nf1
    hold off
    % analysis
    fprintf('When var is '+miui+'\n'+...
    'mean is '+string(mean(data1))+', '+string(mean(data2))+';\n'+...
    'variance is '+string(std(data1))+', '+string(std(data2))+'.');
end





