function [  ] = waveform_plot( y, period_cycle )
%waveform_plot - make a couple plots of different wave positions
% for the body position over a half-period
 
global dim gridsz Gamma

j = 1;
figure(2); clf;
C = linspecer(8); 
for ii=round(period_cycle/2)-round(period_cycle/16):-round(period_cycle/16):1
    yy = y(end-ii,1:gridsz*dim+2);
    plot(yy, 'LineWidth', 4, 'Color', C(j,:)); hold on
    j = j+1;
end; hold off
set(gcf,'color','white')
set(gca, 'LineWidth', 1.0, 'Visible','off');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title(strcat('\Gamma = ',num2str(Gamma)));

end

