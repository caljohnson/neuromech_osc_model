function [ period_cycle, avg_phasediff, wavelength, freq ] = compute_phase_lags( Kappa, ploton)
%compute_phase_lags Computes phase lags between each neuromechanical oscillator
%                   in the chain
%  Uses cross-correlation to get phase lags between each oscillator unit
%  Input: Kappa - gridsz*dim vector of curvatures at each point in body
%          ploton - default 0 for off, pass 1 to plot
%
%  Output: period_cycle - period of oscillation
%          avg_phasediff - avg phase diff of each oscillator to osc. #1
%          wavelength - wavelength of traveling wave
%          freq - frequency of undulation

global Pa TF max_step dim

%if no ploton command passed, default to 0
if nargin < 2
  ploton = 0;
end

%average curvature of each oscillator - Pa is averaging matrix
Kappa = Pa*Kappa(1:end,(TF/max_step)/2:TF/max_step-1);
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
j=1;
for ii=2:dim
   %compute cross-correlation between oscillator ii and 1 
   cross_cor = xcorr(Kappa(ii,:), Kappa(1,:));
   %get index of max-peak of cross-corr
   [~,I_new] = max(cross_cor);
   %phase difference is peak-peak distance / cycle period
   phase_diff(ii) = (I_new-I_control)/period_cycle;
   
   %shift phase differences so that negative phase differences are phase
   %advances, keep track of winding number
%    if phase_diff(ii) <0
%        phase_diff(ii) = phase_diff(ii) + 1;
%    end
   if phase_diff(ii) < phase_diff(ii-1)
       phase_diff(ii) = phase_diff(ii) + j;
   end
   if phase_diff(ii) > j %update winding number
       j = j+1;
   end
end
phase_diff
if ploton==1
    %plot phase differences
    figure(3); clf;
    plot(phase_diff, 'o','MarkerSize',12);
    xlabel('oscillator #');  ylabel('Phase Rel. to Oscillator #1');
end

%compute average phase difference
avg_phasediff = mean(abs(phase_diff(2:end) - phase_diff(1:end-1)));
%compute wavelength of traveling wave
wavelength = 1/(avg_phasediff*dim);
%compute frequency of undulation
freq = 1/(period_cycle*max_step);

end

