function [ init_cond ] = make_init_cond_2pt_case( cycle, T_i, phase_diff)
%make_init_cond Creates the initial condition vector
%   for the pair of oscillators
%   Input: cycle - timeseries of single-oscillator cycle
%          T_i  - period (index-wise) of cycle
%          phase_diff - phase-diff of the pair
%   Output: Vector of initial conditions for each state variable


%if no phase_diff passed, default to 0
if nargin < 3
  phase_diff = 0;
end

%turn 0= < phase_diff <= 1 into an indexed phase-diff
ind_pd = round(T_i*phase_diff)+1;
if ind_pd > T_i %make sure we stay in index array bounds
    ind_pd = T_i;
end

%IC - phase_diff between the two modules
K_0 = zeros(2,1);
AV_0 = zeros(2,1);
AD_0 = zeros(2,1);
volt_V_0  = zeros(2,1);
volt_D_0 = zeros(2,1);

%phase of oscillator 1
    K_0(1) = cycle(1,1);
    AV_0(1) = cycle(1,4);
    AD_0(1) = cycle(1,5);
    volt_V_0(1) = cycle(1,2);
    volt_D_0(1) = cycle(1,3);

%phase of oscillator 2
    K_0(2) = cycle(ind_pd,1);
    AV_0(2) = cycle(ind_pd,4);
    AD_0(2) = cycle(ind_pd,5);
    volt_V_0(2) = cycle(ind_pd,2);
    volt_D_0(2) = cycle(ind_pd,3);

%initial condition
init_cond = [K_0; AV_0; AD_0; volt_V_0; volt_D_0;];

end

