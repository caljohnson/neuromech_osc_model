clear
global t_f t_n t_m c I dt

tees_m = 1:10;

for ll = 1:size(tees_m,2)
% -- NM MODEL PARAMETERS --
t_f=5; t_n = 0.1; t_m = tees_m(ll); %timescales for length, neural, and muscule activity
c = 1; I = 0; %total feedback strength, AVB input bias
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature

% -- coupling parameters --
% proprioceptive coupling strength
eps_prop = 0.1;

%external viscosity (mechanical coupling)
gamma = 1;

% time step size
dt0=0.001;



% ----  I. FIND PERIODIC ORBIT  ----

nmax=100000;  % maximal number of time steps in one period
Kmark=0.1;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt=dt0;
x0(1:nv)=[0;1;-1;0;0]; % initial state (K,V1,V2,A1,A2)

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

% %display result
% figure(1); clf;
%  xlim([0, ii*dt]);
% subplot(4,1,1); plot(v(:,1),'Linewidth', 4); ylabel('\kappa');
% subplot(4,1,2); plot(v(:,2), 'g','Linewidth', 4); hold on; plot(v(:,3), 'r','Linewidth', 4); 
% ylabel('V'); legend('V', 'D');
% subplot(4,1,3); plot(v(:,4), 'g','Linewidth', 4); hold on; plot(v(:,5), 'r','Linewidth', 4); 
% ylabel('A'); legend('V', 'D');
% subplot(4,1,4); plot(tanh(v(:,4)-2)+1, 'g','Linewidth', 4); hold on; plot(tanh(v(:,5)-2)+1, 'r','Linewidth', 4); 
% ylabel('\sigma(A)'); legend('V', 'D');
% suptitle('Cycle timetraces');
% OUTPUT:   v(1:ii,1:nv) contains periodic orbit
%           (ii-1)*dt=~period of orbit


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

% %display result
% figure(2); clf; xlim([0, ii*dt]); suptitle('PRCs');
% subplot(3,1,1); plot(z(:,1), 'Linewidth', 4); ylabel('\kappa');
% subplot(3,1,2); plot(z(:,2), 'g', 'Linewidth', 4); hold on; plot(z(:,3), 'r','Linewidth', 4); 
% ylabel('V'); legend('V', 'D');
% subplot(3,1,3); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
% ylabel('A'); legend('V', 'D');

figure(ll); clf;
xlim([0, ii*dt]);
subplot(2,1,1); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
ylabel('A'); legend('V', 'D'); title(['A_V, A_D PRCs, \tau_m = ' num2str(t_m)]);
subplot(2,1,2); plot(tanh(v(:,4)-2)+1, 'g','Linewidth', 4); hold on; plot(tanh(v(:,5)-2)+1, 'r','Linewidth', 4); 
ylabel('\sigma(A)'); legend('V', 'D'); title('\sigma(A_V), \sigma(A_D) Cycle Timetrace');
end