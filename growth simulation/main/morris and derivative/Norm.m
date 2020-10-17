function [rate,miu] = Norm(Max,Min,miu_num,times,Miu,Sigma)
% 给定最大值、最小值、μ的取样个数和每个μ正态分布的取样的个数
    if Miu==-1&&Sigma==-1
        if miu_num>1
            miu=Min:(Max-Min)/(miu_num-1):Max; % 设定μ的数组
            sigma=zeros(1,miu_num); % 设定σ的数组，每一个μ对应一个σ
            rate=zeros(miu_num,times); % 设定返回的速率的矩阵，一共有miu_num个μ，每一个μ对应times个速率
        else
            miu=Min;
        end
        for i=1:miu_num
            sigma(1,i)=(Max-miu(i))/3; % 认为μ较大时μ+3σ=max。
            if miu-3*sigma<0
                sigma(1,i)=miu/3; % 认为μ较小时μ-3σ=0。
            end
            for j=1:times
                Rate=normrnd(miu(i),sigma(1,i)); % 进行正态分布取值
                if Rate<0 
                    Rate=0; % 小于0的速率纠正为0
                end
%                 for k=1:inf % 小于0的速率重新取值
%                     Rate=normrnd(miu(i),sigma(1,i));
%                     if Rate>0
%                         break;
%                     end
%                 end
                rate(i,j)=Rate; % 对rate赋值
            end
        end
    end
end

