% ==================================================================
%
%                      NM_chainvzg.m
%                      ------- 
%  This program calculates periodic orbits (v), infinitesmal
%  phase response curves (z), and G-functions (g) for the neuromechanical 
%  oscillator chain
%
%   - uses forward Euler method (for now).
%   - calculate the iPRC by linearizing the system around
%     the periodic orbit and solving the adjoint equations.
%   - periodic orbits and iPRCs match those computed with XPP. 
%
%   file "NMchain_step.m" contains FS model equations.
%   file "NMchain_adj.m" contains adjoint equations for linearized 
%         neuromechanical oscillator module model equations.
%   file "gunNMchain.m" calculates G-function.
%
%  CJ 9/30/18
% ==================================================================


clear
global t_f t_n t_m a_v c I dt dim A gamma W

% -- NM MODEL PARAMETERS --
t_f=5; t_n = 0.1; t_m = 1; %timescales for length, neural, and muscule activity
a_v = 1/3; c = 10; I = 0; %voltage scale factor, total feedback strength, AVB input bias
nv=5; % number of variables in oscillator module - 2 neurons, 2 muscles, 1 curvature
dim =6; %number of oscillators in chain

% -- coupling parameters --
% proprioceptive coupling strength
eps_prop = 0.25;
W = speye(dim) - sparse(2:dim, 1:dim-1, eps_prop*ones(dim-1,1), dim,dim);


%external viscosity (mechanical coupling)
gamma = 1e0;
gridsz = 1;
delX = 1;%1/(gridsz*dim+1);
e = ones(gridsz*dim,1);
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);

% %spread operator
% rowinds = 1:gridsz*dim;
% colinds = repelem(1:dim, gridsz);
% S = sparse(rowinds, colinds, ones(gridsz*dim,1));
% 
% %proprioceptive averaging matrix
% rowinds = repelem(1:dim,gridsz);
% colinds = 1:gridsz*dim;
% Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));


% time step size
dt0=0.001;



% ----  I. FIND PERIODIC ORBIT  ----

nmax=100000;  % maximal number of time steps in one period
Kmark=0;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt=dt0;
%initial condition - antiphase alternating
V_ap = repmat([1.7746;-1.7618;], dim/2,1);
A_ap = repmat([1.8033;-1.8134;], dim/2,1);
x0=[0.*ones(dim,1);V_ap; -V_ap; A_ap; -A_ap;]; % initial state (K,V1,V2,A1,A2)

dx=10000.0;     
dxcrit=1e-2 % precision for periodic orbit

tic
while (dx>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
    x=x0;
    y=x;
    ii=1;             
    flagg=0;
    while (flagg==0 && ii<=nmax) % time-step until system returns to
                                 % v=vmark with dv/dt>0
        y=NMchain_step(x,dt); 
        y = y';
        if (x(1)<Kmark && y(1)>=Kmark && ii>1)  
            dt=(Kmark-x(1))/(y(1)-x(1))*dt;
            y=NMchain_step(x,dt);  
            y = y';
            flagg=1; 
            dt=dt0;
        end;
        x=y;
        v(ii,:)=x; %output
        ii=ii+1;
    end;    
    dx=norm(x-x0)
    x0=x;
end;
toc
% rearrange output so that minimum kappa is at t=0.
ii=ii-1;
[vmin,imin]=min(v(:,1));
v2(1:ii-imin+1,:)=v(imin:ii,:);
v2(ii-imin+2:ii,:)=v(1:imin-1,:);
v=v2;

%display result
figure(1); clf;
 xlim([0, ii*dt]);
subplot(4,1,1); plot(v(:,1),'Linewidth', 4); ylabel('\kappa'); hold on
for jj = 2:dim
    plot(v(:,jj),'Linewidth', 4);
end; hold off
subplot(4,1,2); plot(v(:,dim+1),'Linewidth', 4); ylabel('V_V'); hold on
for jj = 2:dim
    plot(v(:,dim+jj),'Linewidth', 4);
end; hold off
subplot(4,1,3); plot(v(:,2*dim+1),'Linewidth', 4); ylabel('V_D'); hold on
for jj = 2:dim
    plot(v(:,2*dim+jj),'Linewidth', 4);
end; hold off
subplot(4,1,4); plot(v(:,3*dim+1),'Linewidth', 4); ylabel('A_V'); hold on
for jj = 2:dim
    plot(v(:,3*dim+jj),'Linewidth', 4);
end; hold off
% subplot(4,1,2); plot(v(:,2), 'g','Linewidth', 4); hold on; plot(v(:,3), 'r','Linewidth', 4); 
% ylabel('V'); legend('V', 'D');
% subplot(4,1,3); plot(v(:,4), 'g','Linewidth', 4); hold on; plot(v(:,5), 'r','Linewidth', 4); 
% ylabel('A'); legend('V', 'D');
% subplot(4,1,4); plot(tanh(v(:,4)-2)+1, 'g','Linewidth', 4); hold on; plot(tanh(v(:,5)-2)+1, 'r','Linewidth', 4); 
% ylabel('\sigma(A)'); legend('V', 'D');
% suptitle('Cycle timetraces');
% OUTPUT:   v(1:ii,1:nv) contains periodic orbit
%           (ii-1)*dt=~period of orbit

return
% ----  II.  CALCULATE iPRC ---- 
%  by solving adjoint of linearized problem (backwards)

y0= ones(nv,1); %initial condition 
dy = 10000;

while (dy>dxcrit) 
% iterate until periodic orbit found to desired precision (dxcrit).
    kk=ii;
    y = y0;
    z(kk,1:nv)=y;            
    while (kk>1)
        kk=kk-1;      
        y=NMadj(y,v(kk,1:nv));
        z(kk,1:nv)=y;
    end;   
    dy=norm(y-y0)
    y0=y;
end;

% normalize Z, i.e. so that dv/dt*z=1 for all t. 
dv=diff(v(1:ii,:))/dt;
z0=(z(1:ii-1,:)+z(2:ii,:))/2;
for k=1:ii-1
sc(k)=z0(k,1:nv)*transpose(dv(k,1:nv));
end;
sc0=median(sc);
z=z/sc0;

%display result
figure(2); clf; xlim([0, ii*dt]); suptitle('PRCs');
subplot(3,1,1); plot(z(:,1), 'Linewidth', 4); ylabel('\kappa');
subplot(3,1,2); plot(z(:,2), 'g', 'Linewidth', 4); hold on; plot(z(:,3), 'r','Linewidth', 4); 
ylabel('V'); legend('V', 'D');
subplot(3,1,3); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
ylabel('A'); legend('V', 'D');


% OUTPUT:   z(1:ii,1:nv) contains iPRC
