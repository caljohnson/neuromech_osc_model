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
TF = 1e4; 
max_step = 1e-2; %time interpolation stepsize
t0 = 0:max_step:TF; %interpolation time steps
Gamma = 1e-1;
c_MAs = linspace(1.3,10,20); 

%mechanical PDE (gamma/kb I + mu/kb AA')k_t = -AA'(k + cM)
%Gamma = gamma/kb, tauf = mu/kb
RHS_matrix = (Gamma*ey+tau_f*(A*A'))\(-(A*A'));

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

%%%%initial phase_diffs
init_phase_diffs = linspace(0,1,21);

%matrix to capture asympototic phase differences
p_diffs = zeros(size(c_MAs,2), size(init_phase_diffs,2)-1);

% ODE simulation code
for ccc = 1:size(c_MAs,2)

    %set feedback strength according to loop
    c_MA = c_MAs(ccc);
    c = c_MA;
    
    %set ODEs
    if thresholding_on == 1
        K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(tanh(AV-2)-tanh(AD-2)));
    else
        K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(AV-AD));
    end

     %compile all RHS ode functions into one        
    ode_rhss = @(t,X) [K_dot(t,X(1:2),X(3:4),X(5:6)); ...
        AV_dot(t,X(7:8),X(9:10),X(3:4));...
         AD_dot(t,X(7:8), X(9:10), X(5:6));...
        volt_V_dot(t, X(1:2), X(7:8));...
        volt_D_dot(t, X(1:2), X(9:10));];


    % ----  I. GET PERIODIC ORBIT  ----
    cycle = vs{ccc};
    ii = size(cycle,1);

    % ODE Simulation Code
    display('running fully coupled ODE simulation');

    %run for diff init phase diffs between the pair
    for mm = 1:size(init_phase_diffs,2)-1
        tic
        %make init condition straight from computed LC
        init_phase_diff = init_phase_diffs(mm);
        init_cond = make_init_cond_2pt_case(cycle, ii, init_phase_diff);
        
        %run simulation
        [t,y] = ode23s(ode_rhss,[0,TF], init_cond);
        %sample cycle at even intervals
        y = interp1(t,y,t0);
        t = t0;
    
        %compute asymptotic phase difference
        Kappa = y(:, 1:2)';
        p_diffs(ccc,mm) = compute_phase_lag_2pt_case(Kappa);
    
        toc
       
    end
end

save('sim_loop_2pt_case_gamma1emin1.mat','p_diffs', 'c_MAs', 'init_phase_diffs');
