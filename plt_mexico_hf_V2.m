% plot mexico 2017 BP results, rupture direction, scale and speed for 3d/1d
% with hf&lf
%
% Author: C. Song, 2017.11.7
%
% Revised based on v1
%

clear;
close all;
figure

% 1d result
subplot (2, 2, 1)
load ('G:\BackProjection\mexico\AL\mexico_al210stations10s0.5HzTo2Hz\al1dhf.mat');
% plot(yuse, xuse, 'go', 'MarkerFaceColor','g', 'markersize', 2); hold on
xs1 = xuse1; 
ys1 = yuse1;
ys2 = yuse2;
xs2 = xuse2;
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
legend([ls1, ls2],'Stage 1','Stage 2', 'location','southwest', 'fontsize', 12);
text(-115, 150, '(a)', 'fontsize', 18, 'color', 'k');
arrow([30, 100], [30, 140], 10, 'Color', [190/255 190/255 190/255], 'Width', 1, 'BaseAngle', 60); hold on
text(26, 146, 'N', 'fontsize', 12, 'color', [190/255 190/255 190/255]);


subplot (2, 2, 3)
vxseq = 0: 0.1: 60;
v1yseq = 2.5*vxseq;
v2yseq = 5*vxseq;
plot(vxseq, v1yseq, 'linestyle', '-', 'color', [190/255 190/255 190/255], 'linewidth', 1); hold on
plot(vxseq, v2yseq, 'linestyle', '-', 'color', [190/255 190/255 190/255], 'linewidth', 1); hold on
text(47, 110, '2.5 km/s', 'fontsize', 11, 'color', [190/255 190/255 190/255], 'fontangle', 'italic');
text(17, 150, '5.0 km/s', 'fontsize', 11, 'color', [190/255 190/255 190/255], 'fontangle', 'italic');
plot(tuse1, distuse1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
plot(tinter11, polyval(coef13, tinter11), 'linestyle', '-', 'color', 'g', 'linewidth', 2); hold on
% xlen = 10;
% arrow([20, 20], [20+xlen, 20+speed13*xlen], 12,'Color','g','Width', 1, 'BaseAngle', 60); hold on
plot(tuse2, distuse2+max(distuse1), 'ko', 'MarkerFaceColor', [132/255 112/255 255/255], 'markersize', 6); hold on
plot(tinter21, polyval(coef23, tinter21)+max(distuse1), 'linestyle', '-', 'color', [132/255 112/255 255/255], 'linewidth', 2); hold on
% dim = [.2 .5 .3 .3];
% annotation('textbox', dim, 'String', 'Stage 1','FitBoxToText','on', 'Color', 'g');
% xlen = 7;
% arrow([45, 20], [45+xlen, 20+speed23*xlen], 12,'Color', [132/255 112/255 255/255], 'Width', 1, 'BaseAngle', 60); hold on
line([36 36], [0 180], 'linewidth', 1, 'color', 'k', 'linestyle', '--'); hold on
text(20, 15, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(45, 15, 'Stage 2', 'fontsize', 12, 'color', 'k');
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 20: 60);
set(gca, 'xlim', [0, 60]);
set(gca, 'ylim', [0, 180]);
set(gca, 'ytick', 0: 20: 180);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('Distance  (km)', 'Fontsize', 15);
text(2, 170, '(c)', 'fontsize', 18, 'color', 'k');



% 3d result
clear;
subplot (2, 2, 2)
load ('G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\al3dhf.mat');
% plot(yuse, xuse, 'go', 'MarkerFaceColor','g', 'markersize', 2); hold on
xs1 = xuse1;
ys1 = yuse1; 
ys2 = yuse2; 
xs2 = xuse2; 
plot(ys1, xs1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
plot(ys2, xs2, 'ko', 'MarkerFaceColor',[132/255 112/255 255/255], 'markersize', 6); hold on
yseq1 = -45: 0.1: 0;             % lat1, lon1 used to plot
xseq1 = polyval(coef11, yseq1);
ls1 = plot(yseq1, xseq1, 'linestyle', '-', 'color', 'g', 'linewidth', 2); hold on
xseq2 = 70: 0.1: 150;             % lat1, lon1 used to plot
yseq2 = polyval(coef21, xseq2);
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
legend([ls1, ls2],'Stage 1','Stage 2', 'location','southwest', 'fontsize', 12);
text(-115, 150, '(b)', 'fontsize', 18, 'color', 'k');
arrow([30, 100], [30, 140], 10, 'Color', [190/255 190/255 190/255], 'Width', 1, 'BaseAngle', 60); hold on
text(26, 146, 'N', 'fontsize', 12, 'color', [190/255 190/255 190/255]);


subplot (2, 2, 4)
vxseq = 0: 0.1: 60;
v1yseq = 2.5*vxseq;
v2yseq = 5*vxseq;
plot(vxseq, v1yseq, 'linestyle', '-', 'color', [190/255 190/255 190/255], 'linewidth', 1); hold on
plot(vxseq, v2yseq, 'linestyle', '-', 'color', [190/255 190/255 190/255], 'linewidth', 1); hold on
text(47, 110, '2.5 km/s', 'fontsize', 11, 'color', [190/255 190/255 190/255], 'fontangle', 'italic');
text(17, 150, '5.0 km/s', 'fontsize', 11, 'color', [190/255 190/255 190/255], 'fontangle', 'italic');
plot(tuse1, distuse1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6); hold on
plot(tinter11, polyval(coef13, tinter11), 'linestyle', '-', 'color', 'g', 'linewidth', 2); hold on
% xlen = 9;
% arrow([20, 20], [20+xlen, 20+speed13*xlen], 12,'Color','g','Width', 1, 'BaseAngle', 60); hold on
plot(tuse2, distuse2+max(distuse1), 'ko', 'MarkerFaceColor', [132/255 112/255 255/255], 'markersize', 6); hold on
plot(tinter21, polyval(coef23, tinter21)+max(distuse1), 'linestyle', '-', 'color', [132/255 112/255 255/255], 'linewidth', 2); hold on
% xlen = 6;
% arrow([50, 20], [50+xlen, 20+speed23*xlen], 12,'Color', [132/255 112/255 255/255], 'Width', 1, 'BaseAngle', 60); hold on
line([37 37], [0 180], 'linewidth', 1, 'color', 'k', 'linestyle', '--'); hold on
text(20, 15, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(45, 15, 'Stage 2', 'fontsize', 12, 'color', 'k');
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 20: 60);
set(gca, 'xlim', [0, 60]);
set(gca, 'ylim', [0, 180]);
set(gca, 'ytick', 0: 20: 180);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('Distance  (km)', 'Fontsize', 15);
text(2, 170, '(d)', 'fontsize', 18, 'color', 'k');


set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 24], 'PaperPositionMode', 'auto');
print('-dpsc','-r600','G:\BackProjection\mexico\AL\mexico_hf.ps');
print('-dpdf','-r600','G:\BackProjection\mexico\AL\mexico_hf.pdf');
