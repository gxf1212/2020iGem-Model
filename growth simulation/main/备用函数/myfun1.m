function dz=myfun1(t,z)
global a21 a12 rho
x = z(1);
y = z(2); 
dz(1) = x*(1-x+a21*y); 
dz(2) = rho*y*(1-y+a12*x); 
dz = [dz(1);dz(2)]; 