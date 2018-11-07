function [ init_cond ] = make_init_cond( cycle, T_i, inphase)
%make_init_cond Creates the initial condition vector
%   for the chain of oscillators
%   Input: cycle - timeseries of single-oscillator cycle
%           
%          inphase - default=1 - in-phase init cond, =0 for antiphase
%   Output: Vector of initial conditions for each state variable

global dim gridsz

%if no inphase passed, default to 1
if nargin < 3
  inphase = 1;
end

if inphase==1
    %IC - each module in-phase
    K_0 = zeros(gridsz*dim,1);
    AV_0 = zeros(dim,1);
    AD_0 = zeros(dim,1);
    volt_V_0  = zeros(dim,1);
    volt_D_0 = zeros(dim,1);
    for ii = 1:dim
            AV_0(ii) = cycle(1,2);
            AD_0(ii) = cycle(1,3);
            volt_V_0(ii) = cycle(1,4);
            volt_D_0(ii) = cycle(1,5);
    end
else
    %IC - antiphase alternating between modules (internally)
    K_0 = zeros(gridsz*dim,1);
    AV_0 = zeros(dim,1);
    AD_0 = zeros(dim,1);
    volt_V_0  = zeros(dim,1);
    volt_D_0 = zeros(dim,1);
    for ii = 1:dim
        if mod(ii,2)==0
            AV_0(ii) = cycle(1,2);
            AD_0(ii) = cycle(1,3);
            volt_V_0(ii) = cycle(1,4);
            volt_D_0(ii) = cycle(1,5);
        else
            AV_0(ii) = cycle(round(T_i/2),2);
            AD_0(ii) = cycle(round(T_i/2),3);
            volt_V_0(ii) = cycle(round(T_i/2),4);
            volt_D_0(ii) = cycle(round(T_i/2),5);
        end
    end
end

%initial condition
init_cond = [K_0; AV_0; AD_0; volt_V_0; volt_D_0;];

end

