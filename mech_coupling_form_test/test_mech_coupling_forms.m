%%%% Test Mech Coupling Forms
% find discrepancy between mechanical coupling forms, if any

global dim gridsz beta cee eps TF Gamma tau_f c_MA tau_m ...
    delX A tau_n a I max_step t0 thresholding_on ey1 ey2

addpath('../waveform_data/');

%simulation runtime
TF = 500;

%global parameters
isold=0; %=1 for old W, =0 for new W
isupstream = 1; %=1 for upstream coupling, =0 for downstream
inphase = 1;%1=start in-phase init cond, 0= antiphase
thresholding_on = 0;
%chain dimension  - number of oscillators
dim = 6;
%improved mechanics factor - gridsize
gridsz = 10;
%proprioception parameters
beta = 0; %ratio of left-right asymmetry in proprioception
cee = 1; %strength of local proprioception
eps = 0; %strength of nonlocal proprioception
%neural parameters
tau_n = 0.001;
a = 1;
I = 0;
%mechanical params
gamma = 1e3; %mech coupling strength via external viscosity/internal visc
kb = 1;
Gamma = gamma/kb;
mu = 5;
tau_f = mu/kb; %mu/kb - internal mechanical timescale - internal viscosity/stiffness
c_MA = 3; %muscle activity mechanical feedback strength = 1/(2adelX) * neural feedback strengths
%timescale of muscle activation
tau_m = 1;
% % Mechanically coupled system 
%make 2nd order centered difference operator A 
%discrete beam PDE matrix with free ends
delX = 1/(gridsz*dim+1);
e = ones(gridsz*dim,1);
% A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], dim, dim+2);
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
%identity matrices
ey1 = speye(gridsz*dim+2,gridsz*dim+2);
ey2 = speye(gridsz*dim,gridsz*dim);

tic
%make ODEs
ode_rhs_ys = mech_coupling_form_odes(1);
ode_rhs_ks = mech_coupling_form_odes(0);

%make init conditions
[ T, cycle, tees ] = get_period_single_oscillator();
T_i = size(tees,2); %get index-length of cycle period

init_cond_ks = make_init_cond(cycle, T_i, inphase);
init_cond_ys = [A\init_cond_ks(1:gridsz*dim);  init_cond_ks(gridsz*dim+1:end);];

%ODE solve stuff
max_step = 1e-1;
t0 = 0:max_step:TF;

[y_ks, ~] = run_model(ode_rhs_ks, init_cond_ks);
[y_ys, t] = run_model(ode_rhs_ys, init_cond_ys);

toc

%%%PLOT KYMOGRAPHS
%plot k form coupling
Kappa = y_ks(:,1:gridsz*dim)';
% FREQUENCY PLOT - plot K vs time
figure(1); clf
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlabel('module #', 'FontSize',30); ylabel('time', 'FontSize',30);
figtitle = ['\epsilon= ', num2str(eps), ...
    ', \beta= ', num2str(beta), ', \Gamma = ', num2str(Gamma,'%.1e')];
title(figtitle, 'FontSize',30);
xlim([1 gridsz*dim+1]); ylim([round(TF/max_step/2) round(TF/max_step)]);
set(gca, 'Xtick', gridsz*(0:dim-1)+gridsz/2);
set(gca, 'XtickLabel', strsplit(num2str(1:dim)));
set(gcf,'color','white')
set(gca, 'LineWidth', 4.0);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(gca, 'FontSize', 30); 
c = colorbar; c.Label.String = 'Curvature';
colormap(fireice);
caxis([-3 3]);
c.Limits = [-1 1];

%plot yform coupling
Kappa = A*y_ys(:,1:gridsz*dim+2)';
% FREQUENCY PLOT - plot K vs time
figure(2); clf
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlabel('module #', 'FontSize',30); ylabel('time', 'FontSize',30);
figtitle = ['\epsilon= ', num2str(eps), ...
    ', \beta= ', num2str(beta), ', \Gamma = ', num2str(Gamma,'%.1e')];
title(figtitle, 'FontSize',30);
xlim([1 gridsz*dim+1]); ylim([round(TF/max_step/2) round(TF/max_step)]);
set(gca, 'Xtick', gridsz*(0:dim-1)+gridsz/2);
set(gca, 'XtickLabel', strsplit(num2str(1:dim)));
set(gcf,'color','white')
set(gca, 'LineWidth', 4.0);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(gca, 'FontSize', 30); 
c = colorbar; c.Label.String = 'Curvature';
colormap(fireice);
caxis([-3 3]);
c.Limits = [-1 1];





