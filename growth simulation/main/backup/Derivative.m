function [derivative_1,derivative_2,Biomass_1,Biomass_2,biomass_1,biomass_2,delta_1,delta_2] = Derivative(Delta_miu,Rate,name)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    global Pf Nf1 Nf2 n1 n2 rc rn rp t_max step
    
    if strcmp(name,"Pf")==1
        Pf=Rate;
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        biomass_1=N1(length(N1)); % 储存最终的生物量
        biomass_2=N2(length(N2));
        Pf=Pf+Delta_miu;
    end
    
    if strcmp(name,"Nf1")==1
        Nf1=Rate;
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        biomass_1=N1(length(N1)); % 储存最终的生物量
        biomass_2=N2(length(N2));
        Nf1=Nf1+Delta_miu;
    end
    
    if strcmp(name,"Nf2")==1
        Nf2=Rate;
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        biomass_1=N1(length(N1)); % 储存第五天的生物量
        biomass_2=N2(length(N2));
        Nf2=Nf2+Delta_miu;
    end
    
    [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
    
    Biomass_1=N1(length(N1)); % 储存第五天的生物量
    Biomass_2=N2(length(N2));
        
    delta_1=Biomass_1-biomass_1; % 储存生物量的变化量
    delta_2=Biomass_2-biomass_2;
    derivative_1=delta_1/Delta_miu;
    derivative_2=delta_2/Delta_miu;     
end

