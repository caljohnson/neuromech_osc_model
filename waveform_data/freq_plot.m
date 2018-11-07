function [ ] = freq_plot( Kappa,saveplot )
%freq_plot Plots curvature kymographs to illustrate freq. of undulation

global dim gridsz beta eps TF Gamma max_step;

%if no saveplot command passed, default to 0
if nargin < 2
  saveplot = 0;
end

% FREQUENCY PLOT - plot K vs time
figure(1); clf
Kappa(gridsz*dim+1,:) = Kappa(gridsz*dim,:);
surf(Kappa');
view(2); shading flat;
xlabel('module #', 'FontSize',30); ylabel('time', 'FontSize',30);
figtitle = ['\epsilon= ', num2str(eps), ...
    ', \beta= ', num2str(beta), ', \Gamma = ', num2str(Gamma,'%.1e')];
title(figtitle, 'FontSize',30);
xlim([1 gridsz*dim+1]); ylim([round(TF/max_step/2) round(TF/max_step)]);
set(gca, 'Xtick', gridsz*(0:dim-1)+gridsz/2);
set(gca, 'XtickLabel', strsplit(num2str(1:dim)));
set(gcf,'color','white')
set(gca, 'LineWidth', 4.0);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(gca, 'FontSize', 30); 
c = colorbar; c.Label.String = 'Curvature';
colormap(fireice);
caxis([-5 5]);
c.Limits = [-3 3];
%save plot
if saveplot==1
    freqplot_title = strcat('freqplot/Gamma_',num2str(Gamma,'%.1e'),...
        'alpha_', num2str(alpha), 'beta', num2str(beta), '.png');
    saveas(gcf, freqplot_title);
end


end

