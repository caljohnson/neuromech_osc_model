
% subroutine to calculate periodic orbit 
%        for model in fsstep
% --------------------------------------


% ----  I. FIND PERIODIC ORBIT  ----

nmax=100000;  % maximal number of time steps in one period
vmark=-75.0;  % initial and final membrane potential (Vm) in periodic orbit (arbitrary!)

dt=dt0;
x0(1:nv)=[-80;0.0;0.5;0.0;0.7]; % initial state (Vm,m,h,n3,n1)

dx=10000.0;     
dxcrit=0.000001 % precision for periodic orbit

while (dx>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
    x=x0;
    y=x;
    ii=1;             
    flagg=0;
    while (flagg==0 && ii<=nmax) % time-step until system returns to
                                 % v=vmark with dv/dt>0
        y=fsstep(x,cbias);      
        if (x(1)<vmark && y(1)>=vmark && ii>1)  
            dt=(vmark-x(1))/(y(1)-x(1))*dt;
            y=fsstep(x,cbias);    
            flagg=1; 
            dt=dt0;
        end;
        x=y;
        v(ii,1:nv)=x(1:nv); %output
        ii=ii+1;
    end;    
    dx=norm(x-x0)
    x0=x;
end;

% rearrange output so that minimum v is at t=0.
ii=ii-1;
[vmin,imin]=min(v(:,1));
v2(1:ii-imin+1,1:nv)=v(imin:ii,1:nv);
v2(ii-imin+2:ii,1:nv)=v(1:imin-1,1:nv);
v=v2;

% OUTPUT:   v(1:ii,1:nv) contains periodic orbit
%           (ii-1)*dt=~period of orbit

