function f=MM(r, k, n)
if length(r)>1 % 为了计算随时间变化的
    f=r./(r+n*k);
else
    f=r/(r+n*k);
end