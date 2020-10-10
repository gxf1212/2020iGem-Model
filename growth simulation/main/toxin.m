function f=toxin(ga, n1, n2)
if length(n1>1)
    f=1-ga*n1.*n2;
else
    f=1-ga*n1*n2;
end