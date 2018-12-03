% plot comparison results after relocation
%
% Last modified: 2018.1.25
% Author: C. Song


clear;
close all;

%% relocation rst
figure;
load('relocation_check.mat');
subplot(2,2,1);
plot(odistuse, orelaplag, 'k*', 'Markersize', 9);
set(gca, 'fontsize', 12);
set(gca, 'xlim', [50, 210]);
set(gca, 'xtick', 50: 50: 205);
set(gca, 'ylim', [-0.1, 0.1]);
set(gca, 'ytick', -0.1: 0.05: 0.1);
xlabel('Distance/(km) ', 'Fontsize', 15);
ylabel('Relative residual', 'Fontsize', 15);
text(55, 0.09, '(a)', 'fontsize', 24, 'color', 'k');

subplot(2,2,2);
plot(ndist, nrelaplag, 'k*', 'Markersize', 9);
set(gca, 'fontsize', 12);
set(gca, 'xlim', [50, 205]);
set(gca, 'xtick', 50: 50: 200);
set(gca, 'ylim', [-0.1, 0.1]);
set(gca, 'ytick', -0.1: 0.05: 0.1);
xlabel('Distance/(km) ', 'Fontsize', 15);
ylabel('Relative residual', 'Fontsize', 15);
text(55, 0.09, '(b)', 'fontsize', 24, 'color', 'k');

%% pws rst
clear;
load('z_shift_stack_reloc.mat');
subplot(2,2,3);
imagesc(tuse1(:, 1), vel, pwssum1'); hold on
colormap(jet);
% colorbar;
caxis([-4.5e-4, 4.5e-4]);
plot(51.2, 8.04, 'k.', 'MarkerSize', 20); hold on
plot(57.09, 8.33, 'k.', 'MarkerSize', 20); hold on
plot(58.96, 5.96, 'k*', 'MarkerSize', 10); hold on
% title('z分量加权相位叠加', 'Fontsize', 18);
set(gca, 'fontsize', 12);
set(gca, 'xtick', 46: 2: 61);
set(gca, 'ytick', 5: 1: 10);
xlabel('Tau/(s) ', 'Fontsize', 15);
ylabel('Velocity/(km·s^{-1}) ', 'Fontsize', 15);
% ylabel(colorbar, '归一化叠加振幅 (cm)');
text(50.5, 7.7, 'Pn', 'fontsize', 18, 'color', 'k');
text(56, 8.0, 'pPn', 'fontsize', 18, 'color', 'k');
text(46.5, 5.2, '(c)', 'fontsize', 24, 'color', 'k');
% set(gcf, 'unit', 'centimeters', 'position', [10, 5, 8, 8]);

clear;
load('t_shift_stack_reloc.mat');
subplot(2,2,4);
imagesc(tuse1(:, 1), vel, pwssum1'); hold on
colormap(jet);
% colorbar;
caxis([-2.5e-4, 2.5e-4]);
plot(86.29, 4.45, 'k.', 'MarkerSize', 20); hold on
plot(96.06, 3.44, 'k*', 'MarkerSize', 10); hold on
% title('t分量加权相位叠加', 'Fontsize', 18);
set(gca, 'fontsize', 12);
set(gca, 'xtick', 76: 4: 98);
set(gca, 'ytick', 2.5: 0.5: 6);
xlabel('Tau/(s) ', 'Fontsize', 15);
ylabel('Velocity/(km·s^{-1}) ', 'Fontsize', 15);
% ylabel(colorbar, '归一化叠加振幅 (cm)');
text(85, 4.2, 'Sn', 'fontsize', 18, 'color', 'k');
text(76.8, 2.65, '(d)', 'fontsize', 24, 'color', 'k');

orient landscape;
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 28, 24], 'PaperPositionMode', 'auto');
print('-dpsc', '-r600', 'relocation_check_rst.ps');
% print('-depsc', '-r600', 'pws_rst.eps');
print('-dpdf', '-r600', 'relocation_check_rst.pdf');
