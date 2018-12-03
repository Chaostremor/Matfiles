% plot mexico 2017 BP results, rupture direction, scale and speed for 3d/1d
% with hf&lf
%
% Author: C. Song, 2017.10.22


clear;
close all;
figure

% 1d result
subplot (2, 2, 1)
load ('G:\BackProjection\mexico\AL\mexico_al210stations10s0.5HzTo2Hz\al1dhf.mat');
% plot(yuse, xuse, 'go', 'MarkerFaceColor','g', 'markersize', 2); hold on
xs1 = [xuse1; x(corind1+1: corind2-1)];
ys1 = [yuse1; y(corind1+1: corind2-1)];
ys2 = [yuse2; y(end2+1: end)];
xs2 = [xuse2; x(end2+1: end)];
plot(ys1, xs1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
plot(ys2, xs2, 'ko', 'MarkerFaceColor', [132/255 112/255 255/255], 'markersize', 6); hold on
yseq1 = -35: 0.1: 0;             % lat1, lon1 used to plot
xseq1 = polyval(coef11, yseq1);
ls1 = plot(yseq1, xseq1, 'linestyle', '-', 'color', 'g', 'linewidth', 2); hold on
yseq2 = -31: 0.1: -25.5;             % lat1, lon1 used to plot
xseq2 = polyval(coef21, yseq2);
ls2 = plot(yseq2, xseq2, 'linestyle', '-', 'color', [132/255 112/255 255/255], 'linewidth', 2); hold on
plot(0, 0, 'rp', 'MarkerFaceColor','r','markersize', 20); hold on
axis equal;
set(gca, 'fontsize', 12);
set(gca, 'XTick', -120: 20: 60);
set(gca, 'xlim', [-120, 60]);
set(gca, 'ylim', [-20, 160]);
set(gca, 'ytick', -20: 20: 160);
xlabel('West-East  (km) ', 'Fontsize', 15);
ylabel('South-North  (km)', 'Fontsize', 15);
legend([ls1, ls2],'Stage 1','Stage 2', 'location','northeast', 'fontsize', 12);
text(-115, 150, '(a)', 'fontsize', 18, 'color', 'k');



subplot (2, 2, 3)
plot(tuse1, distuse1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
xlen = 10;
arrow([20, 20], [20+xlen, 20+speed13*xlen], 12,'Color','g','Width', 1, 'BaseAngle', 60); hold on
plot(tuse2, distuse2, 'ko', 'MarkerFaceColor', [132/255 112/255 255/255], 'markersize', 6); hold on
% dim = [.2 .5 .3 .3];
% annotation('textbox', dim, 'String', 'Stage 1','FitBoxToText','on', 'Color', 'g');
xlen = 5;
arrow([45, 20], [45+xlen, 20+speed23*xlen], 12,'Color', [132/255 112/255 255/255], 'Width', 1, 'BaseAngle', 60); hold on
line([36 36], [0 100], 'linewidth', 1, 'color', 'k', 'linestyle', '--'); hold on
text(20, 10, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(45, 10, 'Stage 2', 'fontsize', 12, 'color', 'k');
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 20: 60);
set(gca, 'xlim', [0, 60]);
set(gca, 'ylim', [0, 100]);
set(gca, 'ytick', 0: 20: 100);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('Distance  (km)', 'Fontsize', 15);
text(2, 95, '(b)', 'fontsize', 18, 'color', 'k');


% 3d result
clear;
subplot (2, 2, 2)
load ('G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\al3dhf.mat');
% plot(yuse, xuse, 'go', 'MarkerFaceColor','g', 'markersize', 2); hold on
xs1 = [xuse1; x(corind1+1: corind2-1)];
ys1 = [yuse1; y(corind1+1: corind2-1)];
ys2 = [yuse2; y(end2+1: end)];
xs2 = [xuse2; x(end2+1: end)];
plot(ys1, xs1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
plot(ys2, xs2, 'ko', 'MarkerFaceColor',[132/255 112/255 255/255], 'markersize', 6); hold on
yseq1 = -45: 0.1: 0;             % lat1, lon1 used to plot
xseq1 = polyval(coef11, yseq1);
ls1 = plot(yseq1, xseq1, 'linestyle', '-', 'color', 'g', 'linewidth', 2); hold on
yseq2 = -53: 0.1: -35;             % lat1, lon1 used to plot
xseq2 = polyval(coef21, yseq2);
ls2 = plot(yseq2, xseq2, 'linestyle', '-', 'color', [132/255 112/255 255/255], 'linewidth', 2); hold on
plot(0, 0, 'rp', 'MarkerFaceColor','r','markersize', 20); hold on
axis equal;
set(gca, 'fontsize', 12);
set(gca, 'XTick', -120: 20: 60);
set(gca, 'xlim', [-120, 60]);
set(gca, 'ylim', [-20, 160]);
set(gca, 'ytick', -20: 20: 160);
xlabel('West-East  (km) ', 'Fontsize', 15);
ylabel('South-North  (km)', 'Fontsize', 15);
legend([ls1, ls2],'Stage 1','Stage 2', 'location','northeast', 'fontsize', 12);
text(-115, 150, '(c)', 'fontsize', 18, 'color', 'k');


subplot (2, 2, 4)
plot(tuse1, distuse1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
xlen = 8;
arrow([20, 20], [20+xlen, 20+speed13*xlen], 12,'Color','g','Width', 1, 'BaseAngle', 60); hold on
plot(tuse2, distuse2, 'ko', 'MarkerFaceColor', [132/255 112/255 255/255], 'markersize', 6); hold on
xlen = 5;
arrow([50, 20], [50+xlen, 20+speed23*xlen], 12,'Color', [132/255 112/255 255/255], 'Width', 1, 'BaseAngle', 60); hold on
line([37 37], [0 100], 'linewidth', 1, 'color', 'k', 'linestyle', '--'); hold on
text(20, 10, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(45, 10, 'Stage 2', 'fontsize', 12, 'color', 'k');
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 20: 60);
set(gca, 'xlim', [0, 60]);
set(gca, 'ylim', [0, 100]);
set(gca, 'ytick', 0: 20: 100);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('Distance  (km)', 'Fontsize', 15);
text(2, 95, '(d)', 'fontsize', 18, 'color', 'k');


set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 24], 'PaperPositionMode', 'auto');
print('-dpsc','-r600','G:\BackProjection\mexico\AL\mexico_hf.ps');
print('-dpdf','-r600','G:\BackProjection\mexico\AL\mexico_hf.pdf');
