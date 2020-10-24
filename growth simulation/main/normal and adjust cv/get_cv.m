function [da1,da2]=get_cv(k,var_name,t_max)

global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global n1 n2 rc rn rp
global max_num_sig times

fprintf('var name=%s \n',var_name);

switch var_name
    case "Pf"
        var=Pf; % record 
    case "Nf1"
        var=Nf1;
    case "r1"
        var=r1;
    case "r2"
        var=r2;
    case "m1"
        var=m1;
    case "m2"
        var=m2;
    case "Kc1"
        var=Kc1;
    case "Kn1"
        var=Kn1;
    case "Kn2"
        var=Kn2;
    case "Kp1"
        var=Kp1;
    case "Kp2"
        var=Kp2;
    case "Q1c"
        var=Q1c;
    case "Q1n"
        var=Q1n;
    case "Q1p"
        var=Q1p;
    case "Q2c"
        var=Q2c;
    case "Q2n"
        var=Q2n;
    case "Q2p"
        var=Q2p;
end

N1_rnd=[];
N2_rnd=[];

% center=20; % sigma=1/center
% sigma=var/(center*2^(k-(max_num_sig+1)/2));

center=20;
sigma=var/center;

re_times=15; % repeat time
d1=zeros(1,re_times);
d2=zeros(1,re_times);

for l=1:re_times
    var_range=normrnd(var,sigma,[1,times]);

    for m=1:times
        switch var_name
            case "Pf"
                Pf=var_range(m);
            case "Nf1"
                Nf1=var_range(m);
            case "r1"
                r1=var_range(m);
            case "r2"
                r2=var_range(m);
            case "m1"
                m1=var_range(m);
            case "m2"
                m2=var_range(m);
            case "Kc1"
                Kc1=var_range(m);
            case "Kn1"
                Kn1=var_range(m);
            case "Kn2"
                Kn2=var_range(m);
            case "Kp1"
                Kp1=var_range(m);
            case "Kp2"
                Kp2=var_range(m);
            case "Q1c"
                Q1c=var_range(m);
            case "Q1n"
                Q1n=var_range(m);
            case "Q1p"
                Q1p=var_range(m);
            case "Q2c"
                Q2c=var_range(m);
            case "Q2n"
                Q2n=var_range(m);
            case "Q2p"
                Q2p=var_range(m);
        end
        % simulation
        [N1, N2, Rc, Rn, Rp, time, G1, G2]=numerical_simulation(n1, n2, rc, rn, rp, t_max);
        N1_rnd=[N1_rnd N1(length(time))];
        N2_rnd=[N2_rnd N2(length(time))];  
    end

    d1(1,l)=std(N1_rnd)/mean(N1_rnd);
    d2(1,l)=std(N2_rnd)/mean(N2_rnd);
end

da1=mean(d1);
da2=mean(d2);

switch var_name
    case "Pf"
        Pf=var;
    case "Nf1"
        Nf1=var;
    case "r1"
        r1=var;
    case "r2"
        r2=var;
    case "m1"
        m1=var;
    case "m2"
        m2=var;
    case "Kc1"
        Kc1=var;
    case "Kn1"
        Kn1=var;
    case "Kn2"
        Kn2=var;
    case "Kp1"
        Kp1=var;
    case "Kp2"
        Kp2=var;
    case "Q1c"
        Q1c=var;
    case "Q1n"
        Q1n=var;
    case "Q1p"
        Q1p=var;
    case "Q2c"
        Q2c=var;
    case "Q2n"
        Q2n=var;
    case "Q2p"
        Q2p=var;
end