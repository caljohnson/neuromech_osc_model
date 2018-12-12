function [ y,t ] = run_model( ode_RHS, init_cond )
%run_model Runs the model simulation
%   for a given ode RHS and init cond

global TF t0;

[t,y] = ode23s(ode_RHS,[0,TF], init_cond);

%sample cycle at even intervals
y = interp1(t,y,t0);
t = t0;


end

