function [ ode_rhss ] = mech_coupling_form_odes(yform)
%MECH_COUPLING_FORM_ODES Creates ODEs for the two mech. coupling forms

global A tau_f Gamma ey1 ey2 c_MA isold isupstream tau_n ...
    tau_m a I thresholding_on gridsz dim

addpath('../waveform_data/');

if yform == 1
    %mechanical PDE (Gamma I + A'A)y_t = -1/tau_f A'(Ay + cM)
    RHS_matrix = (Gamma*ey1+tau_f*(A'*A))\(-A');
    if thresholding_on == 1
        y_dot = @(t,y, AV, AD) RHS_matrix*(A*y + c_MA.*(tanh(AV-2)-tanh(AD-2)));
    else
        y_dot = @(t,y, AV, AD) RHS_matrix*(A*y + c_MA.*(AV-AD));
    end

else
    %mechanical PDE (gamma/kb I + mu/kb AA')k_t = -AA'(k + cM)
    %Gamma = gamma/kb, tauf = mu/kb
    RHS_matrix = (Gamma*ey2+tau_f*(A*A'))\(-(A*A'));
    if thresholding_on == 1
        K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(tanh(AV-2)-tanh(AD-2)));
    else
        K_dot = @(t,K, AV, AD) RHS_matrix*(K + c_MA.*(AV-AD));
    end
end

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
volt_V_dot = @(t,Kappa,volt_V) (1/tau_n).*(volt_V - a.*volt_V.^3 + I + W*Kappa);
volt_D_dot = @(t,Kappa,volt_D) (1/tau_n).*(volt_D - a.*volt_D.^3 + I - W*Kappa);

%muscle eqns: AV' = -AV + (VV-VD)
%             AD' = -AD + (VD-VV)
AV_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_D - volt_V - A);

if yform==1
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
else
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
end


end

