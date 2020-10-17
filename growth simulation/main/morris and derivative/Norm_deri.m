function PF=Norm_deri(miu,times,sigma,mu_step)
PF=zeros(length(miu),times);
% single&same value
if sigma==-1 % evenly spaced delta mu
    time2=round(times/2);
    for j=1:times
        PF(:,j)=abs(j-time2)*mu_step+miu;
    end
end

if mu_step==-1 % normally distributed mu
    for i=1:length(miu)
        PF(i,:)=normrnd(miu(i),sigma,[1,times]);
    end
end

% each miu, each setting
if length(sigma)==length(miu)
    for i=1:length(miu)
            PF(i,:)=normrnd(miu(i),sigma(i),[1,times]);
    end
end

if length(mu_step)==length(miu)
    time2=round(times/2);
    for j=1:times
        PF(:,j)=abs(j-time2)*mu_step+miu;
    end
end
    
end
    