%PRC Loop
% The 2-Point Case study

clear
%Set parameters
global dim gridsz beta cee eps TF Gamma tau_f c_MA tau_m ...
    delX A S Pa tau_n a I max_step t0 thresholding_on ...
    t_f t_n t_m c dt 

addpath('../PRCcode/');
addpath('../waveform_data/');

thresholding_on = 1;
%muscle activity mechanical feedback strength = 1/(2adelX) * neural feedback strengths
c_MAs = linspace(1.3,5,10); 

%chain dimension  - number of oscillators
dim = 2;
%improved mechanics factor - gridsize
gridsz = 1; delX = 1/(gridsz*dim); %grid spacing
%proprioception parameters
beta = 0; %ratio of left-right asymmetry in proprioception
cee = 1; %strength of local proprioception
eps = 0; %strength of nonlocal proprioception
%mechanical params
Gamma = 1e3; %Gamma = gamma/kb - external fluid drag coefficient / stiffness
tau_f = 5; %mu/kb - internal mechanical timescale - internal viscosity/stiffness
%timescale of muscle activation
tau_m = 1;
%continuous voltage model parameters
tau_n = 0.001; a = 1; I = 0;

%constructed globals
%make 2nd order centered difference operator A 
%discrete beam PDE matrix with free ends
e = ones(gridsz*dim,1);
% A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], dim, dim+2);
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%identity
ey = speye(gridsz*dim,gridsz*dim);
%spread operator
rowinds = 1:gridsz*dim;
colinds = repelem(1:dim, gridsz);
S = sparse(rowinds, colinds, ones(gridsz*dim,1));
%proprioceptive averaging matrix
rowinds = repelem(1:dim,gridsz);
colinds = 1:gridsz*dim;
Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));

%redundant, mislabeled params
t_f=tau_f; t_n = tau_n; t_m = tau_m; %timescales for length, neural, and muscule activity
% c = c_MA;

%PRCcode time step size
dt0=1e-4;

%  PRC Code
for ccc = 1:size(c_MAs,2)

    %set feedback strength according to loop
    c_MA = c_MAs(ccc);
    c = c_MA;

dt = dt0; %initially the default time step size
% ----  I. FIND PERIODIC ORBIT  ----

nmax=1e8;  % maximal number of time steps in one period
Kmark=0.2;  % initial and final curvature (K) in periodic orbit (arbitrary!)

dt=dt0;
nv = 5;
x0(1:nv)=[0;1;-1;0;0]; % initial state (K,V1,V2,A1,A2)

dx=10000.0;  
dxcrit=1e-5 % precision for periodic orbit
tic
v = zeros(nmax, nv);
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
toc
% rearrange output so that minimum v is at t=0.
ii=ii-1;
[vmin,imin]=min(v(:,1));
v2(1:ii-imin+1,1:nv)=v(imin:ii,1:nv);
v2(ii-imin+2:ii,1:nv)=v(1:imin-1,1:nv);
v=v2;
clear v2;

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
tic
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
toc
% normalize Z, i.e. so that dv/dt*z=1 for all t. 
dv=diff(v(1:ii,:))/dt;
z0=(z(1:ii-1,:)+z(2:ii,:))/2;
for k=1:ii-1
sc(k)=z0(k,1:nv)*transpose(dv(k,1:nv));
end;
sc0=median(sc);
z=z/sc0;
clear z0 sc;
% %display result
% figure(2); clf; xlim([0, ii*dt]); suptitle('PRCs');
% subplot(3,1,1); plot(z(:,1), 'Linewidth', 4); ylabel('\kappa');
% subplot(3,1,2); plot(z(:,2), 'g', 'Linewidth', 4); hold on; plot(z(:,3), 'r','Linewidth', 4); 
% ylabel('V'); legend('V', 'D');
% subplot(3,1,3); plot(z(:,4), 'g','Linewidth', 4); hold on; plot(z(:,5), 'r','Linewidth', 4); 
% ylabel('A'); legend('V', 'D');

% OUTPUT:   z(1:ii,1:nv) contains iPRC

% ---- III.  CALCULATE G-FUNCTION  ----
% ---------- for the pair --------------
tic
% intrinsic limit cycle info
t=[0:dt:(ii-1)*dt]; 
p=t(ii); % period

% -- Mechanical coupling ---
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
g1=-B(2,1)*(h1m-h2m);
clear h1m h1m;
% [~, Indsg1] = sort(abs(g1));
%prop coupling is not
% g2 =(h1p-h2p);
% [~, Indsh1p] = sort(abs(h1p));

toc
% %%PLOT
% figure(3); clf;
% plot(t(1:ii),g1,'r', 'linewidth',2); hold on; 
% plot(t(1:ii),h1p,'g', 'linewidth',2);
% plot([0,p],[0,0],'k:','linewidth',2);
% plot(t(Indsg1(1:3)),0*Indsg1(1:3),'ko','Markersize',10);
% plot(t(Indsh1p(1:2)),0*Indsh1p(1:2),'ko','Markersize',10);
% legend('mechanical','proprioceptive','"zero"');
% xlabel('time(ms)'); 
% title('G-function for NM paired-oscillator model');

%Save stuff for figures later
vs{ccc} = v;
zs{ccc} = z;
ts{ccc} = t;
g1s{ccc} = g1;
h1ps{ccc} = h1p;
iis{ccc} = ii;

clear v z t g1 h1p ii;

end

save('PRC_loop_2pt_case.mat');
