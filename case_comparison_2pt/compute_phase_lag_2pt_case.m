function [ p_diff ] = compute_phase_lag_2pt_case( Kappa, ploton)
%compute_phase_lags Computes phase lags between the pair

%  Uses cross-correlation to get phase lags between each oscillator unit
%  Input: Kappa - gridsz*dim vector of curvatures at each point in body
%          ploton - default 0 for off, pass 1 to plot
%
%  Output: 
%          p_diff - phase diff of oscillator #2 rel to osc. #1
%          
%          

global TF max_step dim

%if no ploton command passed, default to 0
if nargin < 2
  ploton = 0;
end

%get second-half of time series
Kappa = Kappa(1:end,9*(TF/max_step)/10:TF/max_step-1);
%get autocorrelation of first oscillator
control_corr =  xcorr(Kappa(1,:));
%find peaks of autocorr
[vals,locs] = findpeaks(control_corr);
[~, sorted_locs] = sort(vals,'descend');
%peak-to-peak distance defines period of oscillation
period_cycle = abs(locs(sorted_locs(1))-locs(sorted_locs(2)));
%get max-peak index of autocorr
[~,I_control] = max(control_corr);

%compute phase differences
phase_diff = zeros(dim,1);
winding_no = 0;
for ii=2:dim
   %compute cross-correlation between oscillator ii and 1 
   cross_cor = xcorr(Kappa(ii,:), Kappa(1,:));
   %get index of max-peak of cross-corr
   [~,I_new] = max(cross_cor);
   %phase difference is peak-peak distance / cycle period
   phase_diff(ii) = (I_new-I_control)/period_cycle;
   phase_diff(ii) = phase_diff(ii) + winding_no;  %add winding number
   
   %shift phase differences so that negative phase differences are phase
   %advances, keep track of winding number
   if phase_diff(ii) < phase_diff(ii-1)
       phase_diff(ii) = phase_diff(ii) + 1;
       winding_no = winding_no+1; %update winding number
   end
  
end

if ploton==1
    %plot phase differences
    figure(3); clf;
    plot(phase_diff, 'o','MarkerSize',12);
    xlabel('oscillator #');  ylabel('Phase Rel. to Oscillator #1');
end

%ouput
p_diff = phase_diff(2);

end

