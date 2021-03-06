% calculate cv within normal distributed miu

clear;clc;
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global n1 n2 rc rn rp 

% growth & decay rate
r1=log(2)/30;
r2=0.035/60;
m1=0.196/60; % 14.812 min^{-1} 
m2=0.0035/60;
r2=r2*1.13;

% half-saturation constants, unit: L^{-1}
% B. S
% Kc1=7.4482*0.4; % batch
% Kc1=4.3388*0.4; % fed-batch
Kc1=2.865*0.4;
Kn1=0.638;
Kp1=0.057;
% cyanobacteria
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

% To be realistic, we divide Nostoc into vegetative and N-fixing cells
f_Nf=0.1;f_Ph=1-f_Nf; 

% rates, as constant
% B. S respiration, maintenance energy coefficient of 0.67 mmol glucose/g CDW/h
Re1=0.67e-3*12/60; 
% (pseudomonas aeruginosa) net photosynthesis rate, 4.41+-0.16渭mol(mg路min)
I=100000/54; % luminosity
% NPh=4.41e-3*12.01*tanh(I*0.016/4.41); % 0.0530
NPh=4.41e-3*12.01*I/(165+I+I^2/457); % 0.0103
NPh=NPh*f_Ph;
% bacteria, N decomposition
Nf1=2.600e-4;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/h锛宔stimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
% B.S P-solubilization, inorganic+organic
% Pf=8.87e-4/1440+0; 
Pf=31.0049e-3/1440;


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

rc=rc/1.724; % I forget why?
rn=rn*0.05; % avalible N (inorganic and low-molecular-weight organic)
% rp = rp / 0.06  % just use soluble P

% rho: soil mass / volume in natural state. unit: g/L
% from percentage to mass / volume concentration. unit: g/L
rc=rc*rho;rn=rn*rho;rp=rp*rho; 

% option 2: trial concentration
% rc=0.02;rn=0.001;rp=0.001;rho=1.33e3;


%% adjust
% sample params. you must use "" or it's recognized as chars rather than strings
% var_names=["Pf","Nf1","r1","m1","r2","m2"]; % simple 
var_names=["Pf","Nf1","r1","m1","r2","m2","Q1c","Q1n","Q1p","Q2c","Q2n","Q2p","Kc1","Kn1","Kp1","Kn2","Kp2"];
len=length(var_names);
global max_num_sig times
max_num_sig=1;
days=(1:5)*5;
t_max=1440*days;
times=200;

cv_data=zeros(length(days),len,max_num_sig,2); 

for i=1:length(days)
    for j=1:len
        for k=1:max_num_sig
            fprintf("Program running. i=%d,j=%d,k=%d,",i,j,k);
            [cv_data(i,j,k,1),cv_data(i,j,k,2)]=get_cv(k,var_names(j),t_max(i));
        end
    end
end
% 5*6*5*200*3/1=90000s
% 5*6*1*100*3/1.2=7500s=2h5min
% 5*17*3*200*3/1=153000

%% save data
save('cv_data','cv_data')
% examine
% data5=load('cv_data5%').cv_data;


