% ==================================================================
%
%                      NM_period.m
%                      ------- 
%  This program calculates periodic orbits (v) and computes period
%
%   - uses forward Euler method (for now).
% 
%   file "NMstep.m" contains FS model equations.
%   file "NMadj.m" contains adjoint equations for linearized 
%         neuromechanical oscillator module model equations.
% 
%
%  CJ 11/07/18
% ==================================================================


clear
global t_f t_n t_m c I dt a thresholding_on


cees = 1:0.25:10;
periods = zeros(size(cees,2),1);
max_amp = zeros(size(cees,2),1);

for ll = 1:size(cees,2)

% -- NM MODEL PARAMETERS --
t_f=5; t_n = 0.1; t_m = 1; %timescales for length, neural, and muscule activity
c = cees(ll); I = 0; %total feedback strength, AVB input bias
a = 1;
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature

% -- coupling parameters --
% proprioceptive coupling strength
eps_prop = 0.1;

%external viscosity (mechanical coupling)
gamma = 1e-3;

% time step size
dt0=1e-3;

% ----  I. FIND PERIODIC ORBIT  ----

nmax=1e6;  % maximal number of time steps in one period
Kmark=0.1;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt=dt0;
x0(1:nv)=[0;1;-1;0;0]; % initial state (K,V1,V2,A1,A2)

dx=10000.0;     
dxcrit=1e-5 % precision for periodic orbit

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
periods(ll) = (ii-1)*dt;
max_amp(ll) = max(abs(v(:,1)));
end

% save('NM_periods.mat', 'periods', 'max_amp', 'cees')

figure(1); clf;
subplot(2,1,1); plot(cees, periods, 'o');
xlabel('c'); ylabel('period');
subplot(2,1,2); plot(cees, max_amp, 'o');
xlabel('c'); ylabel('max \kappa');
