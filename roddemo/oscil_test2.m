%
% oscil_test2 -- perform a simulation of an active beam in a viscous
% fluid driven by oscillating momtents in part of the beam.
%
addpath('./src/');  % add path to source directory


Ns = 101;             % number of points
L  = 1;               % length of the beam
ds = L/(Ns-1);        % point spacing
s  = (0:ds:L)';       % arclength coordinate of the beam

ks = 1e4;             % stretching stiffness 
kb = 1e2;              % bending stiffness


Tper = 1;             % period
Nper = 5;            % number of periods to simulate
Tend = Nper*Tper;     % end of simulation
Nt_per_T = 100;       % time steps per period
dt   = Tper/Nt_per_T; % time steps 

mu = 1;               % viscosity
ep = 1.2*ds;          % regularized stokelet parameter

% prescribed curvature function to use
%  the example below is a wave propagating from s=0 to s=L
%
k0  = 5.0;            % curvature amplitude
lam = 2.0;            % wavelength of the wave
kappa_fun = @(s,t)(k0*sin(2*pi/Tper*t)*(s<0.5) );
    
% initialize the swimmer
%
X     = zeros(Ns,2);  % allocate space
th    = zeros(Ns,1);
X(1,1)   = 0;         % initialize so that s=0 is fixed and flat
X(1,2)   = 0;         %   this choice defines the body reference frame
th(1)    = 0;
kappa = kappa_fun(s,0);
for i=2:Ns
    th(i)  = th(i-1)  + ds*kappa(i);
    X(i,1) = X(i-1,1) + ds*cos(th(i-1));
    X(i,2) = X(i-1,2) + ds*sin(th(i-1));
end
X0 = X;


% rotate and recenter -- make the center of mass the origin
%                        make the tip-to-tip angle zero
%
X0 = X0 - repmat(mean(X0),Ns,1);
theta = atan2(X0(Ns,2)-X0(1,2),X0(Ns,1)-X0(1,1));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
X0 = X0*R;

% expand the swimmer to 3D
%
X0 = [X0, zeros(Ns,1)];
  
% make a RHS and Jacibian function for the ODE solver
%
dXdt = @(t,X)velocity3D(X,kappa_fun(s,t),kb,ks,ds,ep,mu);
Jac  = @(t,X)velocity_jac3D(X,kappa_fun(s,t),kb,ks,ds,ep,mu);
opts = odeset('Jacobian',Jac,'RelTol',1e-3,'AbsTol',1e-6);
%opts = odeset('Jacobian',Jac,'RelTol',1e-3,'AbsTol',1e-6,'stats','on');

% call the ode solver
%
tic
[T,Y]=ode15s(dXdt,[0:dt:Tend],X0(:),opts);
toc
% extract the solution from Y
%
Nt = length(T);
Xt = reshape(Y,Nt,Ns,3);
Xt = permute(Xt,[2 3 1]);

% make a movie
%
figure(1);
J = find( s<0.5);
Jp = find( s>=0.5);
Kt = zeros(Ns,Nt);

for k=1:Nt
  plot(Xt(:,1,k),Xt(:,2,k),'b.-',Xt(J,1,k),Xt(J,2,k),'ro');
  axis([-0.7 0.7 -0.5 0.5]);
  pause(0.01);
  
  Kt(:,k) = curvature(Xt(:,1:2,k),ds);
    
end

% plot the kymograph
%
figure(2);
pcolor(T,s,Kt); shading flat;

% plot the mean curvature in the two segments
%
figure(3)
plot(T,mean(Kt(J,:)),T,mean(Kt(Jp,:)));


