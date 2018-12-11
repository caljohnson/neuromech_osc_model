%Simulation Loop:
% The 2-Point Case study

clear
%Set parameters
global dim gridsz beta cee eps TF Gamma tau_f c_MA tau_m ...
    delX A S Pa tau_n a I max_step t0 thresholding_on ...
    t_f t_n t_m c dt 

addpath('../PRCcode/');
addpath('../waveform_data/');

load('PRC_loop_2pt_case.mat');

 %simulation runtime
TF = 2000; 
max_step = 1e-2; %time interpolation stepsize
t0 = 0:max_step:TF; %interpolation time steps

%mechanical PDE (gamma/kb I + mu/kb AA')k_t = -AA'(k + cM)
%Gamma = gamma/kb, tauf = mu/kb
RHS_matrix = (Gamma*ey+tau_f*(A*A'))\(-(A*A'));
if thresholding_on == 1
    K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(tanh(AV-2)-tanh(AD-2)));
else
    K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(AV-AD));
end

%neural coupling matrix
isold=0; %=1 for old W, =0 for new W
isupstream = 1; %=1 for upstream coupling, =0 for downstream
W = make_W_connectivity(isold, isupstream);

volt_V_dot = @(t,Kappa,volt_V) (1/tau_n).*(volt_V - a.*volt_V.^3 + I + W*Kappa);
volt_D_dot = @(t,Kappa,volt_D) (1/tau_n).*(volt_D - a.*volt_D.^3 + I - W*Kappa);

%muscle eqns: AV' = -AV + (VV-VD)
%             AD' = -AD + (VD-VV)
AV_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_D - volt_V - A);

 %compile all RHS ode functions into one        
    % X(1:gridsz*dim) = gridsz*dim K variables
    % X(gridsz*dim+1: gridsz*dim+dim) = dim AV variables
    % X(gridsz*dim+1+dim: gridsz*dim+2*dim) = dim AD variables
    % X(gridsz*dim+1+2*dim: gridsz*dim+3*dim) = dim volt_V variables
    % X(gridsz*dim+1+3*dim: gridsz*dim+4*dim) = dim volt_D variables
    ode_rhss = @(t,X) [K_dot(t,X(1:gridsz*dim),S*X(gridsz*dim+1: gridsz*dim+dim),...
                S*X(gridsz*dim+1+dim: gridsz*dim+2*dim)); ...
        AV_dot(t,X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
        X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1: gridsz*dim+dim));...
         AD_dot(t,X(gridsz*dim+1+2*dim: gridsz*dim+3*dim), ...
        X(gridsz*dim+1+3*dim: gridsz*dim+4*dim),X(gridsz*dim+1+dim: gridsz*dim+2*dim));...
        volt_V_dot(t, Pa*X(1:gridsz*dim), X(gridsz*dim+1+2*dim: gridsz*dim+3*dim));...
        volt_D_dot(t, Pa*X(1:gridsz*dim), X(gridsz*dim+1+3*dim: gridsz*dim+4*dim));];


%%%%initial phase_diffs
init_phase_diffs = linspace(0,1,21);

%matrix to capture asympototic phase differences
p_diffs = zeros(size(c_MAs,2), size(init_phase_diffs,2)-1);

% ODE simulation code
for ccc = 1:size(c_MAs,2)

    %set feedback strength according to loop
    c_MA = c_MAs(ccc);
    c = c_MA;

    % ----  I. GET PERIODIC ORBIT  ----
    cycle = vs{ccc};
    ii = size(cycle,1);

    % ODE Simulation Code
    display('running fully coupled ODE simulation');

   tic
    for mm = 1:size(init_phase_diffs,2)-1
        init_phase_diff = init_phase_diffs(mm);
    %make init condition straight from computed LC
    init_cond = make_init_cond_2pt_case(cycle, ii, init_phase_diff);


    [y,t] = run_model(ode_rhss, init_cond);
   
    Kappa = y(:, 1:gridsz*dim)';
    
    p_diffs(ccc,mm) = compute_phase_lag_2pt_case(Kappa);
    
    end
    toc
end

save('sim_loop_2pt_case.mat','p_diffs', 'c_MAs', 'init_phase_diffs');