% Compare several models
% 2017.9.11
%

clear; close all;
%% model para.
% mod1, Yang Mingzhi, 2007
m1= [4 5.2; 7.5 6.05; 9 6.21; 21 6.55];
mod1= [0 5.2; 4 5.2; 4 6.05;11.5 6.05; 11.5 6.21; 20.5 6.21; 20.5 6.55; 41.5 6.55];

% mod2, crust1.0
lon0 = 107;
lat0 = 38;
b = getcrust(lat0, lon0);
m2(1: 9, 1) = cat(2, b.thk(1: 8), 0);
m2(1: 9, 2) = b.vp(1: 9);
m2(1: 9, 3) = b.vs(1: 9);
m2 = m2(3: end, :);
mod2 = [0 2.5 1.07; 0.5 2.5 1.07; 0.5 4.5 2.51; 4.5 4.5 2.51; 4.5 5 2.88; 7.2 5 2.88; 7.2 6.1 3.55; 23.39 6.1 3.55;...
    23.39 6.3 3.65; 37.74 6.3 3.65; 37.74 7 3.99; 43.99 7 3.99; 43.99 8.25 4.57; 60 8.25 4.57];

% mod4, from inversion, alxa47
m4 = [6.7 5.5 3.28; 31.4 5.96 3.55; 15.6 6.6 4.18; 0 8.12 4.48];
mod4 = [0 5.5 3.28; 6.7 5.5 3.28; 6.7 5.96 3.55; 38.1 5.96 3.55; 38.1 6.6 4.18; 53.7 6.6 4.18; 53.7 8.12 4.48; 60 8.12 4.48];

% mod3, from inversion, alxa45
m3 = [5.8 3.6 2.05; 26.8 6.08 3.64; 16.1 6.5 4.16; 0 8.12 4.48];
mod3 = [0 3.6 2.05; 5.8 3.6 2.05; 5.8 6.08 3.64; 32.6 6.08 3.64; 32.6 6.5 4.16; 48.7 6.5 4.16; 48.7 8.12 4.48; 60 8.12 4.48];

%% calc. travel time curve
h = 24.9;
refdist = 326.0027;
perip = 1.66;
peris = 1.46;
[tpg1, dpg1] = cal_dist_time(4, 4, h, m1(:, 1), m1(:, 2), 'direct');       % mod1, pg

[tpg2, dpg2] = cal_dist_time(7, 5, h, m2(:, 1), m2(:, 2), 'direct');       % mod2, pg
[tpn2, dpn2] = cal_dist_time(7, 5, h, m2(:, 1), m2(:, 2), 'head');       % mod2, pn
[tppn2, dppn2] = cal_dist_time(7, 5, h, m2(:, 1), m2(:, 2), 'ref_head');       % mod2, ppn
[tsg2, dsg2] = cal_dist_time(7, 5, h, m2(:, 1), m2(:, 3), 'direct');       % mod2, sg
[tsn2, dsn2] = cal_dist_time(7, 5, h, m2(:, 1), m2(:, 3), 'head');       % mod2, sn
taupn2 = interp1(dpn2, tpn2, refdist, 'spline')+perip;     % interpolation with linear method
velpn2 = m2(end, 2);
tauppn2 = interp1(dppn2, tppn2, refdist, 'spline')+perip;     % interpolation with linear method
velppn2 = m2(end, 2);
tausn2 = interp1(dsn2, tsn2, refdist, 'spline')+peris;     % interpolation with linear method
velsn2 = m2(end, 3);

[tpg4, dpg4] = cal_dist_time(4, 2, h, m4(:, 1), m4(:, 2), 'direct');       % mod3, pg
[tpn4, dpn4] = cal_dist_time(4, 2, h, m4(:, 1), m4(:, 2), 'head');       % mod3, pn
[tppn4, dppn4] = cal_dist_time(4, 2, h, m4(:, 1), m4(:, 2), 'ref_head');       % mod3, ppn
[tsg4, dsg4] = cal_dist_time(4, 2, h, m4(:, 1), m4(:, 3), 'direct');       % mod3, sg
[tsn4, dsn4] = cal_dist_time(4, 2, h, m4(:, 1), m4(:, 3), 'head');       % mod3, sn
taupn4 = interp1(dpn4, tpn4, refdist, 'spline')+perip;     % interpolation with linear method
velpn4 = m4(end, 2);
tauppn4 = interp1(dppn4, tppn4, refdist, 'spline')+perip;     % interpolation with linear method
velppn4 = m4(end, 2);
tausn4 = interp1(dsn4, tsn4, refdist, 'spline')+peris;     % interpolation with linear method
velsn4 = m4(end, 3);

[tpg3, dpg3] = cal_dist_time(4, 2, h, m3(:, 1), m3(:, 2), 'direct');       % mod3, pg
[tpn3, dpn3] = cal_dist_time(4, 2, h, m3(:, 1), m3(:, 2), 'head');       % mod3, pn
[tppn3, dppn3] = cal_dist_time(4, 2, h, m3(:, 1), m3(:, 2), 'ref_head');       % mod3, ppn
[tsg3, dsg3] = cal_dist_time(4, 2, h, m3(:, 1), m3(:, 3), 'direct');       % mod3, sg
[tsn3, dsn3] = cal_dist_time(4, 2, h, m3(:, 1), m3(:, 3), 'head');       % mod3, sn
taupn3 = interp1(dpn3, tpn3, refdist, 'spline')+perip;     % interpolation with linear method
velpn3 = m3(end, 2);
tauppn3 = interp1(dppn3, tppn3, refdist, 'spline')+perip;     % interpolation with linear method
velppn3 = m3(end, 2);
tausn3 = interp1(dsn3, tsn3, refdist, 'spline')+peris;     % interpolation with linear method
velsn3 = m3(end, 3);

%% plot 

figure
subplot(2, 2, 1);
% vp
plot(mod1(:, 2), mod1(:, 1), 'linewidth', 2, 'color', 'r'); hold on
plot(mod2(:, 2), mod2(:, 1), 'linewidth', 2, 'color', 'g'); hold on
plot(mod3(:, 2), mod3(:, 1), 'linewidth', 2, 'color', [105/255 105/255 105/255]); hold on
plot(mod4(:, 2), mod4(:, 1), 'linewidth', 2, 'color', 'b'); hold on
% vs
plot(mod2(:, 3), mod2(:, 1), 'linewidth', 2, 'color', 'g'); hold on
plot(mod3(:, 3), mod3(:, 1), 'linewidth', 2, 'color', [105/255 105/255 105/255]); hold on
plot(mod4(:, 3), mod4(:, 1), 'linewidth', 2, 'color', 'b'); hold on
line([0 6.55], [41.5 41.5], 'linewidth', 1, 'color', 'r', 'linestyle', '--'); hold on
line([0 7], [43.99 43.99], 'linewidth', 1, 'color', 'g', 'linestyle', '--'); hold on
line([0 6.65], [48.7 48.7], 'linewidth', 1, 'color', [105/255 105/255 105/255], 'linestyle', '--'); hold on
line([0 6.6], [53.7 53.7], 'linewidth', 1, 'color', 'b', 'linestyle', '--'); hold on
h1=legend('杨明芝等(2007)','Crust1.0 ', '初始模型', '最优模型');
% h1=legend('Yang et al., 2007','Crust1.0 ', 'This study');
set(h1, 'location', 'best', 'fontsize', 12);
% set(h1, 'Position', [0.1 0.7 0.2 0.15], 'fontsize', 12);
xlabel('Velocity/(km·s^{-1})', 'Fontsize', 15);
ylabel('Depth/(km) ', 'Fontsize', 15);
set(gca, 'fontsize', 12);
set(gca, 'ydir', 'reverse');       %将z轴反转，使得深度轴竖直向下
set(gca, 'xlim', [0, 9]);
set(gca, 'ylim', [0, 60]);
set(gca, 'xtick', 0: 1: 9);
set(gca, 'ytick', 0: 10: 60);
text(0.3, 4, '(a)', 'fontsize', 24, 'color', 'k');



subplot(2, 2, 2);
% Pg curve
load('alxa45_27_pg_tlag.mat');
plot(pg2, dist, 'k*', 'MarkerSize', 6); hold on;
plot(tpg1, dpg1, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                     % pg mod1
plot(tpg2, dpg2, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;                     % pg mod2
plot(tpg3, dpg3, 'linestyle', '-', 'color', [105/255 105/255 105/255], 'LineWidth', 2); hold on;                     % pg mod3
plot(tpg4, dpg4, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                     % pg mod4
% Sg curve
load('alxa45_27_sg_tlag.mat');
plot(sg2, dist, 'k*', 'MarkerSize', 6); hold on;
plot(tsg2, dsg2, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;                     % pg mod2
plot(tsg3, dsg3, 'linestyle', '-', 'color', [105/255 105/255 105/255], 'LineWidth', 2); hold on;                     % pg mod3
plot(tsg4, dsg4, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                     % pg mod3
h1=legend('直达波到时', '杨明芝等(2007)', 'Crust1.0 ', '初始模型', '最优模型');
% h1=legend('Direct arrival pick', 'Yang et al., 2007', 'Crust1.0 ', 'This study');
% set(h1, 'location', 'best', 'fontsize', 12);
set(h1, 'location', 'southeast', 'fontsize', 12);
xlabel('Travel time/(s) ', 'Fontsize', 15);
ylabel('Distance/(km) ', 'Fontsize', 15);
set(gca, 'fontsize', 12);
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);
text(5, 390, '(b)', 'fontsize', 24, 'color', 'k');



subplot(2, 2, 3);
load('z_shift_stack.mat');
imagesc(tuse1(:, 1), vel, pwssum1'); hold on
colormap(jet);
caxis([-3.5e-4, 3.5e-4]);
plot(taupn2, velpn2, 'kx', 'MarkerSize', 10); hold on
plot(taupn3, velpn3, 'k+', 'MarkerSize', 10); hold on
plot(taupn4, velpn4, 'k*', 'MarkerSize', 10); hold on
plot(tauppn2, velpn2, 'kx', 'MarkerSize', 10); hold on
plot(tauppn3, velppn3, 'k+', 'MarkerSize', 10); hold on
plot(tauppn4, velpn4, 'k*', 'MarkerSize', 10); hold on
set(gca, 'fontsize', 12);
set(gca, 'xtick', 46: 2: 61);
set(gca, 'ytick', 5: 1: 10);
xlabel('Tau/(s) ', 'Fontsize', 15);
ylabel('Velocity/(km·s^{-1})', 'Fontsize', 15);
% h1=legend('crust1.0', '本文');
h1=legend('Crust1.0', '初始模型', '最优模型');
set(h1, 'location', 'northwest', 'fontsize', 12);
text(49, 7.7, 'Pn', 'fontsize', 18, 'color', 'k');
text(55, 7.7, 'pPn', 'fontsize', 18, 'color', 'k');
text(46.5, 9.55, '(c)', 'fontsize', 24, 'color', 'k');


subplot(2, 2, 4);
load('t_shift_stack.mat');
imagesc(tuse1(:, 1), vel, pwssum1'); hold on
colormap(jet);
caxis([-5e-4, 5e-4]);
plot(tausn2, velsn2, 'kx', 'MarkerSize', 10); hold on
plot(tausn3, velsn3, 'k+', 'MarkerSize', 10); hold on
plot(tausn4, velsn4, 'k*', 'MarkerSize', 10); hold on
set(gca, 'fontsize', 12);
set(gca, 'xtick', 76: 4: 98);
set(gca, 'ytick', 2.5: 0.5: 6);
xlabel('Tau/(s) ', 'Fontsize', 15);
ylabel('Velocity/(km·s^{-1})', 'Fontsize', 15);
% h1=legend('crust1.0', '本文');
h1=legend('Crust1.0', '初始模型', '最优模型');
set(h1, 'location', 'northwest', 'fontsize', 12);
% ylabel(colorbar, '归一化叠加振幅 (cm)');
text(85, 4.2, 'Sn', 'fontsize', 18, 'color', 'k');
text(76.8, 5.7, '(d)', 'fontsize', 24, 'color', 'k');

orient landscape;
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 28, 24], 'PaperPositionMode', 'auto');
print('-dpsc','-r600','fig8.ps');
print('-dpdf','-r600','fig8.pdf');


