% plot multi curves in one figure, describing many variables
clear;clc;clf
% parameters are set to global variables
global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p

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
% Nf1=2.65e-4;
% Nf1=0;
% Nostoc N_fixation, 1.72+-0.25 fmol N/cell/h，estimated single cell weight
Nf2=1.72e-15*14.01/60*(7722e6/0.42); 
Nf2=Nf2*f_Nf;
% Nf2=Nf2*70s;
% B.S P-solubilization, inorganic+organic
% Pf=8.87e-4/1440+0; 
% Pf=Pf*10;
Pf=31.0049e-3/1440;
% Pf=Pf*1.5;

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

% call the function
% you should run all cells again to renew parameters
clf;

%% n1,n2
% days=60;
% days=5;
% t_max=1440*days; % 1 day=1440 min
% num=20; % number of lines
% var1=linspace(n1*0.01,n1*2,num); % values of the variable
% var2=linspace(n2*0.01,n2*2,num);

% days=5;
% num=20; % number of lines
% % var=linspace(0.5,5,num); % values of the variable
% % var=linspace(rn,rn,num);
% var=linspace(0.1,2,num);
% 
% Kp=Kp1; % record global variable Kp
% for i=1:num
%     color1=(1/num)*i;
%     color2=(1/num)*(num/1.2-mod(i,round(num/2)));
%     color3=(1/num)*(num-i/3);
%     color=[color1, color2,color3];
%     Kp1=var(i)*Kp;
%     [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
%     
%     figure(1)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, N1,'--', 'linewidth',1.5,'color',color);
%     else 
%         plot(time, N1, 'linewidth', 1.5,'color',color);
%     end 
%     title('growth curve of B.S when days='+string(days))
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     str{i}=("rn="+string(Kp1)); % always need to change!
%     if i==num
% %         legend(str) % matlab is weak in legend typesetting!
%         columnlegend(2,str);
%     end
%     
%     figure(2)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, N2,'--', 'linewidth', 1.5,'color',color);
%     else
%         plot(time, N2, 'linewidth', 1.5,'color',color);
%     end
%     title('growth curve of Nostoc when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     if i==num
% %         legend(str)
%         columnlegend(2,str);
%     end
% end

%% r2,m1
% days=60;
% days=5;
% t_max=1440*days; % 1 day=1440 min
% num=7; % number of lines
% % var=linspace(r2*0.90,r2*1.10,num); % values of the variable
% var=linspace(m1*0.90,m1*1.10,num); 
% for i=1:num
% %     r2=var(i);
%     m1=var(i);
%     [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
%     figure(1)
%     subplot(1,2,1)
%     hold on
%     plot(time, N1, 'linewidth', 1.5);
%     title('growth curve of B.S when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
% %     str{i}=("r2="+string(r2));
%     str{i}=("m1="+string(m1));
%     if i==num
%         legend(str)
%     end
%     
%     subplot(1,2,2)
%     hold on
%     plot(time, N2, 'linewidth', 1.5);
%     title('growth curve of Nostoc when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     if i==num
%         legend(str)
%     end
% end

%% K, Q
clf;
% days=60;
days=20;
t_max=1440*days; % 1 day=1440 min
num=7; % number of lines
% var=linspace(0.5,5,num); % values of the variable
var=linspace(0.6,1.5,num);

% K=Kn1; % record global variable
% Q=Q2n;
ra=Nf1;
% ra=Pf;
for i=1:num
%     color1=(1/num)*i;
%     color2=(1/num)*(num/1.2-mod(i,round(num/2)));
%     color3=(1/num)*(num-i/3);
%     color1=(1/num)*i*0.9;
%     color2=(1/num)*i*0.9;
%     color3=(1/num)*i*0.9;
    color1=(1/num)*(num-i)*0.9;
    color2=(1/num)*(num-i)*0.9;
    color3=1;
    color=[color1, color2,color3];
%     Kn1=var(i)*K;
%     Q2n=var(i)*K;
    Nf1=ra*var(i);
%     Pf=ra*var(i);
    [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
    
    figure()
    if num2str(var(i))=='1'
        plot(time, N1,'--', 'linewidth',1.5,'color',color);
    else 
        plot(time, N1, 'linewidth', 1.5,'color',color);
    end 
    hold on
    title('growth curve of B.S when days='+string(days))
    xlabel('time/min');
    ylabel('biomass concentration/(g\cdot L^{-1})');
%     str{i}=("Pf="+string(Pf)); % always need to change!
    str{i}=("Nf1="+string(Nf1));
    if i==num
%         legend(str) % all in a column
        columnlegend(2,str);
    end
    
    figure()
    if num2str(var(i))=='1'
        plot(time, N2,'--', 'linewidth', 1.5,'color',color);
    else
        plot(time, N2, 'linewidth', 1.5,'color',color);
    end
    hold on
    title('growth curve of Nostoc when days='+string(days));
    xlabel('time/min');
    ylabel('biomass concentration/(g\cdot L^{-1})');
    if i==num
%         legend(str)
        columnlegend(2,str);
    end
    
%     figure(3)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, Rp,'--', 'linewidth', 1.5,'color',color);
%     else
%         plot(time, Rp, 'linewidth', 1.5,'color',color);
%     end
%     title('P when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     if i==num
%         %legend(str)
%         columnlegend(2,str);
%     end
% 
%     figure(4)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, Rn,'--', 'linewidth', 1.5,'color',color);
%     else
%         plot(time, Rn, 'linewidth', 1.5,'color',color);
%     end
%     title('N when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     if i==num
% %         legend(str)
%         columnlegend(2,str);
%     end
    
end

%% resource: the most update and clear figure
% num=20; % number of lines
% % one=linspace(0.9,1,2)(2);
% % var=linspace(0.5,5,num); % values of the variable
% var=linspace(0.1,2,num);
% % var=linspace(rn,rn,num);
% for i=1:num
%     color1=(1/num)*i;
%     color2=(1/num)*(num/1.2-mod(i,round(num/2)));
%     color3=(1/num)*(num-i/3);
%     color=[color1, color2,color3];
%     r0=var(i)*rn;
%     [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, r0, rp, t_max);
%     
%     figure(1)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, N1,'--', 'linewidth',1.5,'color',color);
%     else 
%         plot(time, N1, 'linewidth', 1.5,'color',color);
%     end 
%     title('growth curve of B.S when days='+string(days))
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     str{i}=("rn="+string(r0)); % always need to change!
%     if i==num
% %         legend(str) % matlab is weak in legend typesetting!
%         columnlegend(2,str);
%     end
%     
%     figure(2)
%     hold on
%     if num2str(var(i))=='1'
%         plot(time, N2,'--', 'linewidth', 1.5,'color',color);
%     else
%         plot(time, N2, 'linewidth', 1.5,'color',color);
%     end
%     title('growth curve of Nostoc when days='+string(days));
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
%     if i==num
% %         legend(str)
%         columnlegend(2,str);
%     end
% end

%% P and Nf1
% num=7;
% var=linspace(Nf1*0.95,Nf1*1.05,num);
% for i=1:num
%     Nf1=var(i);
%     [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
%     figure(1)
%     hold on
%     ti=sprintf('Nf1=%e',Nf1);
%     title(['f of B.S when ',ti]);
%     subplot(3,3,i)
%     pl3_3f=plot(time, G1./N1/r1, 'o','MarkerSize', 1.5);
%     pl3_3r=plot(time, MM(Rc, Kc1, N1), time, MM(Rn, Kn1, N1), time, MM(Rp, Kp1, N1));
%     legend('f','C','N','P')
%     xlabel('time/min');
%     figure(2)
%     hold on
%     plot(time, Rp, 'linewidth', 1.5);
%     title('P concentration');
%     str{i}=[ti];
%     xlabel('time/min');
%     ylabel('biomass concentration/(g\cdot L^{-1})');
% end
% legend(str)
