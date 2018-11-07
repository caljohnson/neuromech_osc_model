function [ T, cycle, tees] =get_period_single_oscillator_v2( tau_f, c_MA, tau_m, ploton )
%get_period_single_oscillator_v2
%gets the period and limit cycle of the single, uncoupled oscillator
%with continuous neural dynamics

%if no ploton command passed, default to 0
if nargin < 4
  ploton = 0;
end

%simulation runtime
TF = 1e2;

kappa_dot = @(t,kappa,AV,AD) (1/tau_f)*(-kappa - c_MA*(tanh(AV-2) - tanh(AD-2)));
   
%continuous voltage model
tau_n = 0.001;
a = 1/3;
I = 0;
volt_V_dot = @(t,Kappa,volt_V) (1/tau_n).*(volt_V - a.*volt_V.^3 + I + Kappa);
volt_D_dot = @(t,Kappa,volt_D) (1/tau_n).*(volt_D - a.*volt_D.^3 + I - Kappa);


%muscle eqns: AV' = -AV + (VV-VD)
%             AD' = -AD + (VD-VV)
AV_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/tau_m).*(volt_D - volt_V - A);

%IC 
K_0 = -1;
AV_0 = -2;
AD_0 = 2;
volt_V_0 = 1;
volt_D_0 = 0;

%ODE solve stuff
max_step = 1e-2;
t0 = 0:max_step:TF;

%compile all RHS ode functions into one        
ode_rhss = @(t,X) [kappa_dot(t,X(1),X(2),X(3)); ...
    AV_dot(t,X(4), X(5), X(2));...
    AD_dot(t,X(4), X(5), X(3));...
    volt_V_dot(t, X(1), X(4));...
    volt_D_dot(t, X(1), X(5));];
init_cond = [K_0; AV_0; AD_0; volt_V_0; volt_D_0;];

[t,y] = ode23(ode_rhss,[0,TF], init_cond);
    
%sample cycle at even intervals
y = interp1(t,y,t0);
t = t0';

% plot(y(round((3/4)*size(t,1)):size(t,1),1),'o');
start_flag = 0;
tol = 1e-1;
for ii=round((3/4)*size(t,1)):size(t,1)
   if start_flag==1
       if t(ii) < start_time + 10*max_step;
           continue
       end
   end
   if abs(y(ii,1)) < tol && start_flag == 0 
       start_time = t(ii);
       start_ind = ii;
       start_flag = 1;
       if y(ii,1)>y(ii-1,1)
           upward_flag = 1;
       else
           upward_flag = 0;
       end
   elseif abs(y(ii,1)) < tol && start_flag == 1 && y(ii,1)>y(ii-1,1) && upward_flag == 1
           end_time = t(ii);
           end_ind = ii;
           break
   elseif abs(y(ii,1)) < tol && start_flag == 1 && y(ii,1)<y(ii-1,1) && upward_flag == 0
            end_time = t(ii);
            end_ind = ii;
            break
   end
end

T = end_time - start_time;
tees = t(start_ind:end_ind) - start_time;
return_inds = [1 2 3 4 5];
cycle = y(start_ind:end_ind,return_inds);

if ploton==1
    figure(10); clf;
    subplot(2,3,1); plot(tees, cycle(:,1), '.'); ylabel('\kappa'); xlabel('t');%ylim([-1 1]);
    subplot(2,3,2); plot(tees, cycle(:,2), '-'); ylabel('A_V'); xlabel('t');
    subplot(2,3,3); plot(tees, cycle(:,3), '-'); ylabel('A_D'); xlabel('t');
    subplot(2,3,4); plot(tees, cycle(:,4), '-'); ylabel('V^V'); xlabel('t');
    subplot(2,3,5); plot(tees, cycle(:,5), '-'); ylabel('V^D'); xlabel('t');
end

end

