function [ T, cycle, tees] =get_period_single_oscillator(ploton)
%get_period_single_oscillator
%gets the period and limit cycle of the single, uncoupled oscillator
%with continuous neural dynamics

global dt

addpath('../PRCcode/');

%if no ploton command passed, default to 0
if nargin < 1
  ploton = 0;
end

% ----  I. FIND PERIODIC ORBIT  ----
nv = 5;
nmax=1e8;  % maximal number of time steps in one period
Kmark=0.1;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt0=1e-3;
dt = dt0;
x0(1:nv)=[0;1;-1;0;0]; % initial state (K,V1,V2,A1,A2)

dx=10000.0;     
dxcrit=1e-5; % precision for periodic orbit

while (dx>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
    x=x0;
    y=x;
    ii=1;             
    flagg=0;
    while (flagg==0 && ii<=nmax) % time-step until system returns to
                                 % v=vmark with dv/dt>0
        y=NMstep(x);      
        if (x(1)<Kmark && y(1)>=Kmark && ii>1)  
            dt=(Kmark-x(1))/(y(1)-x(1))*dt;
            y=NMstep(x);    
            flagg=1; 
            dt=dt0;
        end;
        x=y;
        v(ii,1:nv)=x(1:nv); %output
        ii=ii+1;
    end; 
    if ii == nmax
        display('no periodic solution found');
        exit;
    end
    dx=norm(x-x0);
    x0=x;
end;

% rearrange output so that minimum v is at t=0.
ii=ii-1;
[~,imin]=min(v(:,1));
v2(1:ii-imin+1,1:nv)=v(imin:ii,1:nv);
v2(ii-imin+2:ii,1:nv)=v(1:imin-1,1:nv);
v=v2;

cycle = v;
T = (ii-1)*dt;
tees = 0:dt:T;

if ploton==1
    figure(10); clf;
    subplot(2,3,1); plot(tees, v(:,1), '.'); ylabel('\kappa'); xlabel('t');%ylim([-1 1]);
    subplot(2,3,2); plot(tees, v(:,4), '-'); ylabel('A_V'); xlabel('t');
    subplot(2,3,3); plot(tees, v(:,5), '-'); ylabel('A_D'); xlabel('t');
    subplot(2,3,4); plot(tees, v(:,2), '-'); ylabel('V^V'); xlabel('t');
    subplot(2,3,5); plot(tees, v(:,3), '-'); ylabel('V^D'); xlabel('t');
end
end

