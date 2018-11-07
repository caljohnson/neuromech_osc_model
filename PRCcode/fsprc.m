

% subroutine to calculate iPRC for the model in fsstep
% given the periodic orbit v calculated by fsper
% ----------------------------------------------------


% ----  II.  CALCULATE iPRC ---- 
%  by solving adjoint of linearized problem (backwards)

y(1:nv)=[-0.007,-0.003,5.0,-6.0,2.0]; % initial condition 

for j=1:10  % iterate jmax=10 times for convergence to Z. 
    kk=ii;
    z(kk,1:nv)=y(1:nv);
    while (kk>1)
        kk=kk-1;      
        y=fsadj(y,v(kk,1:nv));
        z(kk,1:nv)=y;
    end;
end;

% normalize Z, i.e. so that dv/dt*z=1 for all t. 
dv=diff(v(1:ii,:))/dt;
z0=(z(1:ii-1,:)+z(2:ii,:))/2;
for k=1:ii-1
sc(k)=z0(k,1:nv)*transpose(dv(k,1:nv));
end;
sc0=median(sc);
z=z/sc0;

% OUTPUT:   z(1:ii,1:nv) contains iPRC


