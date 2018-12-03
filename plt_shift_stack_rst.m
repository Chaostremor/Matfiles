% plot free shift stack
%
% Last modified: 2018.1.7
% Author: C. Song

% %% pws rst
% clear;
% figure;
% load('z_shift_stack.mat');
% subplot(1,2,1);
% imagesc(tuse1(:, 1), vel, pwssum1'); hold on
% colormap(jet);
% % colorbar;
% caxis([-3.5e-4, 3.5e-4]);
% plot(50.54, 8.12, 'k.', 'MarkerSize', 20); hold on
% plot(56.48, 7.98, 'k.', 'MarkerSize', 20); hold on
% plot(58.73, 6, 'k*', 'MarkerSize', 10); hold on
% % title('z分量加权相位叠加', 'Fontsize', 18);
% set(gca, 'fontsize', 12);
% set(gca, 'xtick', 46: 2: 61);
% set(gca, 'ytick', 5: 1: 10);
% xlabel('Tau/(s) ', 'Fontsize', 15);
% ylabel('Velocity/(km·s^{-1}) ', 'Fontsize', 15);
% % ylabel(colorbar, '归一化叠加振幅 (cm)');
% text(49, 7.8, 'Pn', 'fontsize', 18, 'color', 'k');
% text(55, 7.7, 'pPn', 'fontsize', 18, 'color', 'k');
% text(46.5, 5.2, '(a)', 'fontsize', 24, 'color', 'k');
% % set(gcf, 'unit', 'centimeters', 'position', [10, 5, 8, 8]);
% 
% clear;
% load('t_shift_stack.mat');
% subplot(1,2,2);
% imagesc(tuse1(:, 1), vel, pwssum1'); hold on
% colormap(jet);
% % colorbar;
% caxis([-5e-4, 5e-4]);
% plot(85.96, 4.48, 'k.', 'MarkerSize', 20); hold on
% plot(96.39, 3.69, 'k*', 'MarkerSize', 10); hold on
% % title('t分量加权相位叠加', 'Fontsize', 18);
% set(gca, 'fontsize', 12);
% set(gca, 'xtick', 76: 4: 98);
% set(gca, 'ytick', 2.5: 0.5: 6);
% xlabel('Tau/(s) ', 'Fontsize', 15);
% ylabel('Velocity/(km·s^{-1}) ', 'Fontsize', 15);
% % ylabel(colorbar, '归一化叠加振幅 (cm)');
% text(85, 4.3, 'Sn', 'fontsize', 18, 'color', 'k');
% text(76.8, 2.65, '(b)', 'fontsize', 24, 'color', 'k');
% 
% orient landscape;
% set(gcf, 'unit', 'centimeters', 'position', [10, 8, 28, 10], 'PaperPositionMode', 'auto');
% print('-dpsc', '-r600', 'pws_rst.ps');
% % print('-depsc', '-r600', 'pws_rst.eps');
% print('-dpdf', '-r600', 'pws_rst.pdf');


%% phase check
figure
clear;
load('z_shift_stack.mat');
subplot(1, 2, 1);
for i = 1: nsum
    temp1 = data2(:, i)./(max(abs(data2(:, i)))) .*30.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t2(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', 'k', 'LineWidth', 0.8); hold on;
    fill(t2(:, i), temp2+ double(dist(i+index-1)), [156/255 156/255 156/255], 'edgealpha', 0); hold on;
end
set(gca, 'fontsize', 12);
set(gca, 'XTick', 30: 10: 80);
set(gca, 'xlim', [30, 80]);
set(gca, 'ylim', [240, 420]);
set(gca, 'ytick', 250: 50: 400);
xlabel('Travel time/(s) ', 'Fontsize', 15);
ylabel('Distance/(km)', 'Fontsize', 15);
distcurve1 = 250: 0.1: 400;
p1 = 1/8.12;         % Pn
tau1 = 50.54;
y1 = p1.*(distcurve1 - refdist)+tau1;
plot(y1, distcurve1, 'linestyle', '-', 'color', 'b', 'LineWidth', 1.5); hold on;
distcurve2 = 340: 0.1: 400;
p2 = 1/7.98;                   % pPn
tau2 = 56.48;
y2 = p2.*(distcurve2 - refdist)+tau2;     
plot(y2, distcurve2, 'linestyle', '-', 'color', 'r', 'LineWidth', 1.5); hold on;
text(58, 410, 'Pn', 'fontsize', 18, 'color', 'b');
text(64, 410, 'pPn', 'fontsize', 18, 'color', 'r');
text(32, 410, '(c)', 'fontsize', 24, 'color', 'k');


load('t_shift_stack.mat');
subplot(1, 2, 2);
for i = 1: nsum
    temp1 = data2(:, i)./(max(abs(data2(:, i)))) .*30.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t2(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', 'k', 'LineWidth', 0.8); hold on;
    fill(t2(:, i), temp2+ double(dist(i+index-1)), [156/255 156/255 156/255], 'edgealpha', 0); hold on;
end
set(gca, 'fontsize', 12);
set(gca, 'XTick', 60: 10: 110);
set(gca, 'xlim', [60, 110]);
set(gca, 'ylim', [240, 420]);
set(gca, 'ytick', 250: 50: 400);
xlabel('Travel time/(s) ', 'Fontsize', 15);
ylabel('Distance/(km)', 'Fontsize', 15);
distcurve1 = 250: 0.1: 400;
p1 = 1/4.48;         % Sn
tau1 = 85.96;
y1 = p1.*(distcurve1 - refdist)+tau1;
plot(y1, distcurve1, 'linestyle', '-', 'color', 'b', 'LineWidth', 1.5); hold on;
text(100, 410, 'Sn', 'fontsize', 18, 'color', 'b');
text(62, 410, '(d)', 'fontsize', 24, 'color', 'k');

orient landscape;
set(gcf, 'unit', 'centimeters', 'position', [10, 8, 28, 10], 'PaperPositionMode', 'auto');
print('-dpsc','-r600','pws_wave_phase.ps');
% print('-depsc','-r600','pws_wave_phase.eps');
print('-dpdf','-r600','pws_wave_phase.pdf');
