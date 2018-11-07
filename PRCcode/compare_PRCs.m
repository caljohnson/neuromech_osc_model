global t_f t_n t_m c I dt a

cees = 2:0.125:3;

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

nmax=1e8;  % maximal number of time steps in one period
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

% %Plot PRC for A and sigma(A)
% figure(ll); clf;
% xlim([0, ii*dt]);
% subplot(2,1,1); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
% ylabel('A'); legend('V', 'D'); title(['A_V, A_D PRCs, c = ' num2str(c)]);
% subplot(2,1,2); plot(tanh(v(:,4)-2)+1, 'g','Linewidth', 4); hold on; plot(tanh(v(:,5)-2)+1, 'r','Linewidth', 4); 
% ylabel('\sigma(A)'); legend('V', 'D'); title('\sigma(A_V), \sigma(A_D) Cycle Timetrace');
%

%Plot PRCs
figure(ll); clf; xlim([0, ii*dt]); suptitle(['PRCs, c = ' num2str(c)]);
subplot(3,1,1); plot(z(:,1), 'Linewidth', 4); ylabel('\kappa');
subplot(3,1,2); plot(z(:,2), 'g', 'Linewidth', 4); hold on; plot(z(:,3), 'r','Linewidth', 4); 
ylabel('V'); legend('V', 'D');
subplot(3,1,3); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
ylabel('A'); legend('V', 'D');
saveas(figure(ll), ['compare_prc_plots/PRCs_c_' num2str(c) '.png']);

% OUTPUT:   z(1:ii,1:nv) contains iPRC

% ---- III.  CALCULATE G-FUNCTION  ----
% ---------- for the pair --------------

% intrinsic limit cycle info
t=[0:dt:(ii-1)*dt]; 
p=t(ii); % period

% -- Mechanical coupling ---
dim = 2; %pair of 2 units
%second-difference matrix A
A = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)], [0 1 2], dim, dim+2);
B = inv(full(A*A'));
%K' from limit cycle
Kp = dv;
Kp(ii+1) = Kp(1);
%mech coupling terms
% e = (A*A')\(gamma*Kp);

% calculate H-functions 
% ---------------------
for j=1:ii
    %mechanical coupling
  m1 = z(j:ii,1).*transpose(Kp(1:ii+1-j));
  m2 = z(1:j-1,1).*transpose(Kp(ii+2-j:ii));
    %h function - mech
  h1m(j)=(sum(m1)+sum(m2))*dt0/p;
  h2m(ii+1-j)=(sum(m1)+sum(m2))*dt0/p;
    %proprioceptive coupling
  p1=z(j:ii,2).*v(1:ii+1-j,1);
  p2=z(1:j-1,2).*v(ii+2-j:ii,1);
    %h function - prop
  h1p(j)=(sum(p1)+sum(p2))*dt0/p;
  h2p(ii+1-j)=(sum(p1)+sum(p2))*dt0/p;
end;

% calculate G-functions for the Pair
% ---------------------
%mech coupling is symmetric in the pair!
g1=B(2,1)*(h1m-h2m);
ming1 = min(abs(g1));
%prop coupling is not
% g2 =(h1p-h2p);
minh1p = min(abs(h1p));
figure(ll); clf;
plot(t(1:ii),g1,'r', 'linewidth',2); hold on; 
plot(t(1:ii),h1p,'g', 'linewidth',2);
plot([0,p],[0,0],'k:','linewidth',2);
plot(t(abs(g1)<=1e-4),0*t(abs(g1)<=1e-4),'ko','Markersize',10);
plot(t(abs(h1p)<=1e-4),0*t(abs(h1p)<=1e-4),'ko','Markersize',10);
legend('mechanical','proprioceptive','"zero"');
xlabel('time(ms)'); 
title(['G-functions for NM paired-oscillator model, c = ' num2str(c)]);
saveas(figure(ll), ['compare_g_fns/g_fns_c' num2str(c) '.png']);



clear v v2 sc z t h1m h2m h1p h2p
end