% main for growth simulation
% real parameters, simple and full simulation
% unit of time is min^{-1} in accord with those of rates
% rates per biomass, Ks are deminsionless?
% we assume that the unit of N is g/L

% MATLAB supports versions Python 2.7, 3.6, and 3.7

clear;clc;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0

% growth & decay rate
r1=log(2)/30;
r2=0.035/60;
m1=0.196/60; % 14.812 min^{-1} 
m2=0.0035/60;
r2=r2*1.13;

% m1=m1*0.95;
% r2=r2*0.34;
% m2=m2*0.34;

% half-saturation constants, unit: L^{-1}
% B. S
% Kc1=7.4482*0.4; % batch
% Kc1=4.3388*0.4; % fed-batch
Kc1=2.865*0.4;
Kn1=0.638;
Kp1=0.057;
% cyanobacteria
% Kn2=9.31;
Kn2=0.050993; % lansijun
Kp2=0.018739; % tonglvjiadanbao

% Q, constant
% B. S
% Q1c=0.303; % estimated elemental mass fraction
Q1c=0.34*8/0.3*1e-3*72; % reciprocal of growth/biomass yield
Q1n=0.0734;
Q1p=0.00607;
% Nostoc
Q2c=0.2523; % estimated elemental mass fraction
%Q2c=0.488/0.65; % not sensitivity, little change in result
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
I=100000/54; % luminosity on a sunny day
% I=45.7; % literature value
% NPh=4.41e-3*12.01*tanh(I*0.016/4.41); % 0.0530
NPh=4.41e-3*12.01*I/(165+I+I^2/457); % 0.0103
NPh=NPh*f_Ph;
% bacteria, N decomposition
Nf1=2.600e-4;
% Nf1=2.5e-4;
% Nf1=2.202e-4;
% Nf1=0;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/h，estimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
% Nf2=Nf2*70s;
% B.S P-solubilization, inorganic+organic
% Pf=8.87e-4/1440+0; 
% Pf=Pf*10;
Pf=31.0049e-3/1440;

% pH
pH0=9.1;


%% simulation setup
% variables
% biomass, 1 for B.S, 2 for Nostoc
% c,n,p for the element

% initial values
% cell density, Nostoc > B.S often
n1=0.1;n2=0.2;
% nutrient concentration, unit: g/L

% option 1: natural condition
% sand dune
rc=7.2e-3;rn=2.6e-4;rp=8.88e-7;rho=1.426e3;
% plants
% rc=2.35e-2;rn=2.8e-4;rp=3.676e-6;rho=1.446e3;

rc=rc/1.724; % organic carbon
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

% eps params
alpha1=1.4071;
beta1=0.0093;
alpha2=10.1682;
beta2=0.0057;
% ratio=2;alpha1=alpha1*ratio;beta1=beta1*ratio; % project enhancement
params=[alpha1 beta1 alpha2 beta2];

% days=60;
days=0.5;
t_max=1440*days; % 1 day=1440 min

% set searching range
% num1=10;num2=10;
% var1=linspace(Nf1*0.9,Nf1*1.00,num1); % Nf1
% var2=linspace(Pf/4,Pf/4*5,num2); % Pf
% var1=linspace(rn*0.6,rn*1.5,num1); % rn
% var2=linspace(rp*1,rp*10,num2);

% var1=linspace(n1*0.1,n1*2,num1); % n1
% var2=linspace(n2*0.1,n2*2,num2);
num1=10;num2=10;
var1=linspace(n1*0.6,n1*1.5,num1); % Nf1
var2=linspace(n2*0.6,n2*1.5,num2);

%% 2d matrix with no-decrease..., v2
% using matrix search area, each row for a sample
data=zeros(num1*num2,5); 
% n1, n2, total eps, no-decrease, resource
for i=1:num1
    for j=1:num2
        fprintf('i=%d, j=%d\n',i,j);
        idx=(i-1)*num1+j;
        n1=var1(i);
        n2=var2(j);
        [N1, N2, Rc, Rn, Rp, time, G1, G2,E1,E2]=numerical_simulation_eps(n1, n2, rc, rn, rp, t_max,params);
        
        data(idx,1)=n1;data(idx,2)=n2;
        data(idx,3)=E1(t_max)+E2(t_max);
        data(idx,4)=and(N1(720)>n1,N2(720)>n2);
    end
end

% processing
idx=find(data(:,4)==1);
data_pro=data(idx,:);
data_pro=norm_col(data_pro);
 
coef=[-1;-1;1;0;1];
evaluate=data_pro*coef;


%% with EPS, v1
% biomass=zeros(num1,num2,2);
% EPS=zeros(num1,num2,2);
% for i=1:num1
%     for j=1:num2
%         fprintf('i=%d, j=%d\n',i,j);
%         n1=var1(i);
%         n2=var2(j);
%         [N1, N2, Rc, Rn, Rp, time, G1, G2,E1,E2]=numerical_simulation_eps(n1, n2, rc, rn, rp, t_max,params);
%         biomass(i,j,1)=N1(t_max);
%         biomass(i,j,2)=N2(t_max);
%         EPS(i,j,1)=E1(t_max);
%         EPS(i,j,2)=E2(t_max);
%     end
% end
% 
% % plot
% clc;clf;
% 
% [x,y]=meshgrid(var1,var2);
% figure(1)
% % BS
% subplot(1,2,1)
% mesh(x,y,biomass(:,:,1));
% hold on
% % scatter3(var1(5),var2(1),biomass(5,1,1),'*');
% title('B.S biomass')
% xlabel('var1')
% ylabel('var2')
% hold on
% % idx=30; % fix n1 and see the curve
% % % n10=var1(idx); 
% % % plot3(x(:,idx),y(:,idx),biomass(:,idx,1),'r*');
% % n20=var2(idx); 
% % plot3(x(idx,:),y(idx,:),biomass(idx,:,1),'r*');
% % zlim([0,0.1]);
% 
% % Nostoc
% subplot(1,2,2)
% mesh(x,y,biomass(:,:,2));
% % hold on
% % plot3(var1(5),var2(1),biomass(5,1,2));
% title('Nostoc biomass')
% xlabel('var1')
% ylabel('var2')
% % zlim([0.6,0.8]);
% 
% figure(2)
% % BS
% subplot(1,2,1)
% mesh(x,y,EPS(:,:,1));
% hold on
% % scatter3(var1(5),var2(1),biomass(5,1,1),'*');
% title('B.S EPS')
% xlabel('var1')
% ylabel('var2')
% 
% % Nostoc
% subplot(1,2,2)
% mesh(x,y,EPS(:,:,2));
% title('Nostoc EPS')
% xlabel('var1')
% ylabel('var2')
% 
% evaluate=zeros(num1,num2);
% coef=[1;1;-2;-2];
% for i=1:num1
%     for j=1:num2
%         evaluate(i,j)=[EPS(i,j,1) EPS(i,j,2) var1(i) var2(j)]*coef;
%     end
% end
% figure(3)
% mesh(x,y,evaluate)

