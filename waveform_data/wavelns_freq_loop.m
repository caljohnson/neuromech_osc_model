%wavelns_freq_loop.m
% BIG PARAMETER LOOP to collect data on waveform and freq of traveling wave
%simulates mechanical coupling on discrete neuron-system chain
% with beam equation model, plus extra neural coupling
% and with "nondimensionalized" parameters

global dim gridsz alpha beta cee TF Gamma tau_f c_MA tau_m ...
    delX A S Pa tau_n a I max_step t0 eps

% W params
isold=0; %=1 for old W, =0 for new W
isupstream=1; %=1 for upstream coupling, =0 for downstream

%init cond param
inphase = [1,0]'; %=1 for inphase init cond, =0 for antiphase

%chain dimension  - number of oscillators
dim = 6;
%improved mechanics factor - gridsize
gridsz = 10;

% Gammas = 1e2;
% taufs = 5;

Gammas = logspace(-4,4,20)';%[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]';
% taufs = 5;%[2 3 4 5 8 10 15 20 50 100]';

wavelengths = zeros(size(Gammas,1), size(inphase,1));

for aa = 1:size(Gammas,1)
    for dd = 1:size(inphase,1)
        
%asymmetry parameters
alpha = 0.1; %ratio of nonlocal to local coupling
beta = 0; %ratio of left-right asymmetry
cee = 1; %strength of total neuromechanical inputs
eps = alpha;

%simulation runtime
TF = 200;
tic;    
%mechanical params
Gamma = Gammas(aa); %mech coupling strength via external viscosity/internal visc
tau_f = 5; %mu/k - internal mechanical timescale - internal viscosity/stiffness
c_MA = 30; %muscle activity mechanical feedback strength = 1/(2adelX) * neural feedback strengths
%timescale of muscle activation
tau_m = 1;

% % Mechanically coupled system 
%make 2nd order centered difference operator A 
%discrete beam PDE matrix with free ends
delX = 1/(gridsz*dim+1);
e = ones(gridsz*dim,1);
% A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], dim, dim+2);
A = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);

%identity
ey = speye(gridsz*dim+2,gridsz*dim+2);

%mechanical PDE (Gamma I + A'A)y_t = -1/tau_f A'(Ay + cM)
RHS_matrix = (Gamma*ey+A'*A)\(-A'/tau_f);
y_dot = @(t,y, AV, AD) RHS_matrix*(A*y + c_MA.*(tanh(AV-2)-tanh(AD-2)));

%spread operator
rowinds = 1:gridsz*dim;
colinds = repelem(1:dim, gridsz);
S = sparse(rowinds, colinds, ones(gridsz*dim,1));

%proprioceptive averaging matrix
rowinds = repelem(1:dim,gridsz);
colinds = 1:gridsz*dim;
Pa = sparse(rowinds, colinds, (1/gridsz).*ones(gridsz*dim,1));

%neural coupling matrix
W = make_W_connectivity(isold, isupstream);

%continuous voltage model
tau_n = 0.001;
a = 1/3;
I = 0;
volt_V_dot = @(t,Kappa,volt_V) (1/tau_n).*(volt_V - a.*volt_V.^3 + I + W*Kappa);
volt_D_dot = @(t,Kappa,volt_D) (1/tau_n).*(volt_D - a.*volt_D.^3 + I - W*Kappa);

%muscle eqns: AV' = -AV + (VV-VD)
%             AD' = -AD + (VD-VV)
AV_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_D - volt_V - A);

%pick initial conditions from cycle
%simulate single-oscillator to get limit cycle and pick initial phases
[ T, cycle, tees ] = get_period_single_oscillator_v2( tau_f, c_MA, tau_m );
T_i = size(tees,1); %get index-length of cycle period


%start inphase
init_cond = make_init_cond(cycle, inphase(dd));

%ODE solve stuff
max_step = 1e-2;
t0 = 0:max_step:TF;

%compile all RHS ode functions into one        
% X(1:gridsz*dim+2) = gridsz*dim+2 y variables - gridsz*dim Kappa =  A*y
% X(gridsz*dim+3: gridsz*dim+2+dim) = dim AV variables
% X(gridsz*dim+3+dim: gridsz*dim+2+2*dim) = dim AD variables
% X(gridsz*dim+3+2*dim: gridsz*dim+2+3*dim) = dim volt_V variables
% X(gridsz*dim+3+3*dim: gridsz*dim+2+4*dim) = dim volt_D variables
ode_rhss = @(t,X) [y_dot(t,X(1:gridsz*dim+2),S*X(gridsz*dim+3: gridsz*dim+2+dim),...
            S*X(gridsz*dim+3+dim: gridsz*dim+2+2*dim)); ...
    AV_dot(t,X(gridsz*dim+3+2*dim: gridsz*dim+2+3*dim), ...
    X(gridsz*dim+3+3*dim: gridsz*dim+2+4*dim),X(gridsz*dim+3: gridsz*dim+2+dim));...
     AD_dot(t,X(gridsz*dim+3+2*dim: gridsz*dim+2+3*dim), ...
    X(gridsz*dim+3+3*dim: gridsz*dim+2+4*dim),X(gridsz*dim+3+dim: gridsz*dim+2+2*dim));...
    volt_V_dot(t, Pa*A*X(1:gridsz*dim+2), X(gridsz*dim+3+2*dim: gridsz*dim+2+3*dim));...
    volt_D_dot(t, Pa*A*X(1:gridsz*dim+2), X(gridsz*dim+3+3*dim: gridsz*dim+2+4*dim));];

tic
[y,t,Kappa] = run_model(ode_rhss, init_cond);
toc;


%compute waveform info and plot
[period_cycle, avg_phasediff, wavelength, freq] = compute_phase_lags(Kappa, 0);

wavelengths(aa,dd) = wavelength;
    end
end

data_title = strcat('asc_newW_wavelns_freq_loop_eps_', num2str(eps),'.mat');
save(data_title);


%plot stuff
% figure(1); clf;
% plot(log10(Gammas), wavelengths, '-o', 'MarkerSize', 12);
% lg_labels1 = {size(taufs,1)};
% for ii=1:size(taufs,1)
%     lg_labels1(ii) = {strcat('\tau_f = ', num2str(taufs(ii)))};
% end
% legend(lg_labels1);
% xlabel('log_{10}\Gamma'); ylabel('Wavelength');

% figure(2); clf;
% plot(taufs, freqs, '-o', 'MarkerSize', 12);
% lg_labels2 = {size(Gammas,1)};
% for ii=1:size(Gammas,1)
%     lg_labels2(ii) = {strcat('\Gamma = ', num2str(Gammas(ii),'%.1e'))};
% end
% legend(lg_labels2);
% xlabel('\tau_f'); ylabel('frequency');



