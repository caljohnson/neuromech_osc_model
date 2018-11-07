% ==================================================================
%
%                      figures_Loop_NMvzg.m
%                      ------- 
%  This program draws the periodic orbits (v) and  infinitesmal
%  phase response curves (z) from the giant parameter loop
%
%
%  CJ 10/01/18%
% ==================================================================

cees = 1:20;
tees_m = 0.1:0.5:5;

periods = zeros(size(cees,2), size(tees_m,2));

for kk = 1:size(cees,2)
    for ll = 1:size(tees_m,2)
         
% -- NM MODEL PARAMETERS --
t_f=5; t_n = 0.1; t_m = tees_m(ll); %timescales for length, neural, and muscule activity
c = cees(kk); I = 0; %total feedback strength, AVB input bias
nv=5; % number of variables in model - 2 neurons, 2 muscles, 1 curvature

load(['loopdata_box/tm_' num2str(t_m,'%.1e') '_c_' num2str(c) '.mat'])

% cs(kk,ll) = c;
% ts_m(kk,ll) = t_m;
periods(kk,ll) = (ii-1)*dt;

    end
end

save('loopdata_box/periods.mat', 'cees', 'tees_m', 'periods');

%surf(tees_m(2:end),cees,periods(:,2:end))
% xlabel('t_m'); ylabel('c');