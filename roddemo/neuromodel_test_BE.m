%
% neuromodel_test_BE -- perform a simulation of an active beam in a viscous
% fluid driven by neuromechanics (no proprioceptive coupling)
%
addpath('./src/');  % add path to source directory
% addpath('../waveform_data/');

Ns = 60;             % number of points
L  = 1;               % length of the beam
ds = L/(Ns-1);        % point spacing
s  = (0:ds:L)';       % arclength coordinate of the beam

ks = 1e4;             % stretching stiffness 
kb = 1e-2;              % bending stiffness


%module stuff
c = 30; %module feedback strength
no_mods = 6; %number of neuromechanical module units
no_mod_points = Ns/no_mods; %number of body points per unit 
%spread operator
rowinds = 1:Ns;
colinds = repelem(1:no_mods, no_mod_points);
S = sparse(rowinds, colinds, ones(Ns,1));
%proprioceptive averaging matrix
rowinds = repelem(1:no_mods,no_mod_points);
colinds = 1:Ns;
Pa = sparse(rowinds, colinds, (1/no_mod_points).*ones(Ns,1));


Tper = 1;             % period
Nper = 100;            % number of periods to simulate
Tend = Nper*Tper;     % end of simulation
Nt_per_T = 100;       % time steps per period
dt   = Tper/Nt_per_T; % time steps 
Nt = Nt_per_T*Nper;

mu = 1;               % viscosity
ep = 1.2*ds;          % regularized stokelet parameter

% % prescribed curvature function to use
% %  the example below is a wave propagating from s=0 to s=L
% %
k0  = 5.0;            % curvature amplitude
% lam = 2.0;            % wavelength of the wave
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
neuro_module = repmat([1 -1 0 0],1,no_mods);


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

T = [0:dt:Tend];
Xt = zeros(Ns,3,Nt+1);
Xt(:,:,1)= X0;
X = X0;
tic
for k =1:Nt
   t = k*dt;
   %insert time step for internal dynamics
%    kappa = kappa_fun(s,t);  %replace with active moment * kb
   neuro_module = NM_step3d(neuro_module, X, Pa, dt, no_mods);
   driving_moment = moment_step(neuro_module,c, no_mods, S);
   X = BE_step3D(X,driving_moment,dt,kb,ks,ds,ep,mu);
   Xt(:,:,k+1)=X;
end
toc

% % make a movie
% 
% figure(1);
% J = find( s<0.5);
% Jp = find( s>=0.5);
% Kt = zeros(Ns,Nt);
% 
% for k=1:Nt
%  plot(Xt(:,1,k),Xt(:,2,k),'b.-',Xt(J,1,k),Xt(J,2,k),'b.-');
%   axis([-0.7 0.7 -0.5 0.5]);
%   pause(0.01);
%   
%  Kt(:,k) = curvature(Xt(:,1:2,k),ds);
%     
% end
% 
% % plot the kymograph
% 
% figure(2);
% pcolor(T(1:end-1),s,Kt); shading flat;
% 
% plot the mean curvature in the two segments
% 
% figure(3)
% plot(T,mean(Kt(J,:)),T,mean(Kt(Jp,:)));

save('test.mat');

