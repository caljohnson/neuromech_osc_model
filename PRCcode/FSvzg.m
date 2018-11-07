
% ==================================================================
%
%                      FSvzg.m
%                      ------- 
%  This program calculates periodic orbits (v), infinitesmal
%  phase response curves (z), and G-functions (g) for the Erisir
%  etal. (1999) model for cortical fast-spiking (FS) interneurons.
%
%   - uses forward Euler method (for now).
%   - calculate the iPRC by linearizing the system around
%     the periodic orbit and solving the adjoint equations.
%   - periodic orbits and iPRCs match those computed with XPP. 
%
%   file "fsstep.m" contains FS model equations.
%   file "fsadj.m" contains adjoint equations for linearized 
%         FS model equations.
%   file "gunFS.m" calculates G-function.
%
%  TL 9/29/03
% ==================================================================

clear
global gNa gK1 gK3 gL VNa VK VL cap dt

% -- FS MODEL PARAMETERS --
gNa=112.5; gK3=225.0; gK1=0.225; gL=0.25; 
VNa=74.0; VK=-90.0; VL=-70.0; cap=1.0;
nv=5; % number of variables in model
cbias=10.0; % applied current

% -- coupling parameters --
% rate and strength of alpha-synapse 
a=0.75; gs=-2.0;
% strength of gap junctional conductance
gc=0.1;

% heterogeneity in I_app(cbias)
dIapp=1.0;


% time step size
dt0=0.001;



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





% ---- III.  CALCULATE G-FUNCTION  ----

% intrinsic limit cycle info
t=[0:dt:(ii-1)*dt]; 
p=t(ii) % period

% alpha-current inhibitory synaptic current 
% ------------------------------------------
% p-periodic and triggered when v=0.0,v'>0. 
s=exp(-a*p);
for i=1:ii
e0(i)=a^2*exp(-a*t(i))*(t(i)*(1.0-s)+p*s)/((1.0-s)^2);
end;
for i=2:ii
if (v(i,1)>0.0 & v(i-1,1)<0.0)
     istar=ii-i;
end;
end;
e(1:ii-istar)=e0(istar+1:ii);
e(ii-istar+1:ii)=e0(1:istar);


% plot v(t), Z(t)=adj(v(t)) and synaptic current
% ----------------------------------------------
figure(1);
clf;
subplot(311);
plot(t(1:ii),v(:,1),'linewidth',2);
hold on;
plot(t(1:ii)+p,v(:,1),'r','linewidth',2);
plot(t(1:ii)-p,v(:,1),'g','linewidth',2);
axis([-p,2*p,-100,70]);
ylabel('v');
subplot(312);
plot(t(1:ii),z(:,1),'linewidth',2);
hold on;
plot(t(1:ii)+p,z(:,1),'r','linewidth',2);
plot(t(1:ii)-p,z(:,1),'g','linewidth',2);
axis([-p,2*p,-0.1,0.4]);
ylabel('v-adj');
subplot(313);
plot(t(1:ii),e(:),'linewidth',2);
hold on;
plot(t(1:ii)+p,e(:),'r','linewidth',2);
plot(t(1:ii)-p,e(:),'g','linewidth',2);
axis([-p,2*p,0,0.5]);
xlabel('time (ms)');
ylabel('I-syn');
pause(1.0);



% calculate H-functions 
% ---------------------
for j=1:ii
r1=z(j:ii,2).*transpose(e(1:ii+1-j));
r2=z(1:j-1,2).*transpose(e(ii+2-j:ii));
h1s(j)=(sum(r1)+sum(r2))*dt0/p;
h2s(ii+1-j)=(sum(r1)+sum(r2))*dt0/p;
r1=z(j:ii,1).*(v(1:ii+1-j,1)-v(j:ii,1));
r2=z(1:j-1,1).*(v(ii+2-j:ii,1)-v(1:j-1,1));
h1c(j)=(sum(r1)+sum(r2))*dt0/p;
h2c(ii+1-j)=(sum(r1)+sum(r2))*dt0/p;
end;

% ... and G-functions
% -------------------
g1=gs*(h1s-h2s);
g2=gc*(h1c-h2c);


% effects of heterogeneity in I_app
sz=dIapp*sum(z(:,1))*dt/p;

% plot G-functions for different rho
% ----------------------------------
% "animated"
for rho=0:0.1:1.0
g=(1-rho)*g1+rho*g2;
figure(2)
clf;
plot(t(1:ii),g(:),'linewidth',2);
hold on;
plot([0,p],[sz,sz],'k:','linewidth',2);
plot(t(1:ii),g2(:),'r-','linewidth',2);
plot(t(1:ii),g1(:),'g-','linewidth',2);
legend('inhibit+electric','"zero"','electric','inhibit');
xlabel('time(ms)'); 
str1=sprintf('G-function for FS model; alpha= %0.4g, rho= %0.4g',a,rho);
title(str1);
pause(1.0);
end;



