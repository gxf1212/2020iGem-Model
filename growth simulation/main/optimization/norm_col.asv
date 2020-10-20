function data_norm=norm_col(data_pro)
% feature normalization
data_norm=data_pro;
siz=size(data_pro);
for i=1:siz(2)
    mu=mean(data_pro(:,i));
    sig=std(data_pro(:,i));
    if and(sig==0,i~=4)
        fprintf("sigma=0 when i=%d. return mu in this column.\n",i);
        continue;
    end
    data_norm(:,i)=(data_norm(:,i)-mu)/sig;
end