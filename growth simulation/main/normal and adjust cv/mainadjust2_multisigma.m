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

% toxinic supression, gamma, unit: (g/L)^{-2}
ga1=1; % 2 on 1
ga2=1;
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
% (pseudomonas aeruginosa) net photosynthesis rate, 4.41+-0.16æ¸­mol(mgè·¯min)
I=100000/54; % luminosity
% NPh=4.41e-3*12.01*tanh(I*0.016/4.41); % 0.0530
NPh=4.41e-3*12.01*I/(165+I+I^2/457); % 0.0103
NPh=NPh*f_Ph;
% bacteria, N decomposition
Nf1=2.600e-4;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/hé”›å®”stimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
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

rc=rc/1.724; % I forget why?
rn=rn*0.05; % avalible N (inorganic and low-molecular-weight organic)
% rp = rp / 0.06  % just use soluble P

% rho: soil mass / volume in natural state. unit: g/L
% from percentage to mass / volume concentration. unit: g/L
rc=rc*rho;rn=rn*rho;rp=rp*rho; 

% option 2: trial concentration
% rc=0.02;rn=0.001;rp=0.001;rho=1.33e3;

%% adjust
% sample params
% the value we pick before
% mu=Pf;
%mu=r1;
%mu=r2;
mu=m1;
%mu=m2;
%mu=2.600e-4; %Nf1
times=100;
days=10;
t_max=1440*days; % 1 day=1440 min
fsave1=[];fsave01=[];
fsave2=[];fsave02=[];
 [N10, N20, Rc0, Rn0, Rp0, time0, G10, G20]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
for i=1:1
%     sigma=mu/(5*2^(i-1));   %i=3,sigma=mu*5%
%     sigma=mu*0.1;
    sigma=mu*0.05;
%     sigma=mu*0.025;
    range=normrnd(mu,sigma,[1,times]);
    t=1:times;
    N1_rnd=[];
    N2_rnd=[];
% simulation params
for j=1:times
    j
    m1=range(j);
    [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
    N1_rnd=[N1_rnd N1(length(time))];
    N2_rnd=[N2_rnd N2(length(time))];  
end
aver1=mean(N1_rnd)    %average value
aver2=mean(N2_rnd) 
std1=std(N1_rnd)    %standard deviation
std2=std(N2_rnd)
cv1=std1/aver1
cv2=std2/aver2     %coefficient of variation
 [h1,p1]=lillietest( N1_rnd)   %whether the data obtained conforms to normal distribution , if h=0 ¡°Yes¡±
 [h2,p2]=lillietest( N2_rnd)
%[H1,sig1,ci1] = ttest(N1_rnd,aver1,0.05,1) 
%[H2,sig2,ci2] = ttest(N2_rnd,aver2,0.05,1)
%sig is the probability of the observed value
%ci is the 1-alpha confidence interval with true mean value. 
%H = 0 means the original hypothesis cannot be rejected at the significance level alpha; 
%% analysis
% if i==3
% figure(1)
% scatter(t,N1_rnd)
% title('Scatter diagram of B.S');
% figure(2)
% scatter(t,N2_rnd)
% title('Scatter diagram of Nostoc');
% figure(3)
% %[a1,b1]=hist(N1_rnd);   
% %bar(b1,a1/sum(a1));  % Frequency histogram
% histfit(N1_rnd)
% title('Frequency histogram of B.S ');
% xlabel('Biomass of B.S');
% ylabel('Frequency');
% figure(4)
% %[a2,b2]=hist(N2_rnd);
% %bar(b2,a2/sum(a2)); 
% histfit(N2_rnd)
% title('Frequency histogram of Nostoc ');
% xlabel('Biomass of Nostoc');
% ylabel('Frequency');
% end
fsave01=[aver1 std1 cv1 h1];
fsave1=[fsave1;fsave01]
fsave02=[aver2 std2 cv2 h2];
fsave2=[fsave2;fsave02]
end
% Ne1= N10(t_max)
% Ne2= N20(t_max)
%save('file1.mat','N1_rnd','aver1','std1','cv1','h1'); 
%save('file2.mat','N2_rnd','aver2','std2','cv2','h2'); 