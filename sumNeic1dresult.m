% plot mexico 2017 BP results, rupture direction, scale and speed for 1d
% NEIC location
%
% Author: C. Song, 2018.11.18
%
%

%%     NEIC LOCATION
clear;
close all;
figure;


% 1d result
clear;
subplot (1, 2, 1)
load ('/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_al_neic210stations10s0.5HzTo2Hz/al1dhf.mat');
% plot(yuse, xuse, 'go', 'MarkerFaceColor','g', 'markersize', 2); hold on
xs1 = xuse1;
ys1 = yuse1; 
ys2 = yuse2; 
xs2 = xuse2;
h1 = plot(0, 0, 'rp', 'MarkerFaceColor','r','markersize', 20); hold on
h2 = plot(ys1, xs1, 'ko', 'MarkerFaceColor', 'g', 'markersize', 6);
h3 = plot(ys2, xs2, 'ko', 'MarkerFaceColor',[132/255 112/255 255/255], 'markersize', 6);
yseq1 = -47: 0.1: 0;             % lat1, lon1 used to plot
xseq1 = polyval(coef11, yseq1);
% xseq1 = 0: 0.1: max(x1)-1;             % lat1, lon1 used to plot
% yseq1 = polyval(coef11, xseq1);
ls1 = plot(yseq1, xseq1, 'linestyle', '-', 'color', 'g', 'linewidth', 2);
% xseq2 = min(x2): 0.1: max(x2)+4;             % lat1, lon1 used to plot
% yseq2 = polyval(coef21, xseq2);
% yseq2 = min(y2): 0.1: y2(1);             % lat1, lon1 used to plot
% xseq2 = polyval(coef21, yseq2);
yseq2 = -38: 0.1: -33;             % lat1, lon1 used to plot
xseq2 = polyval(coef21, yseq2);
ls2 = plot(yseq2, xseq2, 'linestyle', '-', 'color', [132/255 112/255 255/255], 'linewidth', 2);
axis equal;
set(gca, 'fontsize', 12);
set(gca, 'XTick', -120: 20: 60);
set(gca, 'xlim', [-120, 60]);
set(gca, 'ylim', [-20, 160]);
set(gca, 'ytick', -20: 20: 160);
xlabel('West-East  (km) ', 'Fontsize', 15);
ylabel('South-North  (km)', 'Fontsize', 15);
legend([ls1, ls2],'Stage 1','Stage 2','location','southwest', 'fontsize', 12, 'AutoUpdate', 'off');
text(-115, 150, '(a)', 'fontsize', 18, 'color', 'k');
h5 = text(26, 146, 'N', 'fontsize', 12, 'color', [190/255 190/255 190/255]);
set(gca, 'position', [0.1 0.3 0.35 0.35]);
load('/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_al_neic210stations10s0.5HzTo2Hz/parret.mat');
latmean = mean(ret.lat);
lonmean = mean(ret.lon);
[darc, azi] = distance(ret.lat0, ret.lon0, latmean, lonmean);
xmean = deg2km(latmean - ret.lat0);
ymean = deg2km((lonmean - ret.lon0)*cosd(ret.lat0));
angle=atand(xmean/ymean);
if coef11(1)>=0       % rupture to right, look from the lower plane
    azi1 = angle;
else                           %  rupture to left
    azi1 = 360+angle;
end
dy=20;
h7 = text(-110, 95, 'Array', 'fontsize', 12, 'color', [160/255 160/255 160/255]);
h8 = text(-40, 140, '1-D', 'fontsize', 14, 'color', [100/255 100/255 100/255]);
h6 = arrow([-70, 90], [-70-dy, 90+dy*tand(azi1-270)], 12, 'Color', [160/255 160/255 160/255], 'Width', 2, 'BaseAngle', 60);
h4 = arrow([30, 100], [30, 140], 10, 'Color', [190/255 190/255 190/255], 'Width', 1, 'BaseAngle', 60);

subplot (1, 2, 2)
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
line([37.3 37.3], [0 180], 'linewidth', 1, 'color', 'k', 'linestyle', '--'); hold on
text(20, 15, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(45, 15, 'Stage 2', 'fontsize', 12, 'color', 'k');
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 20: 60);
set(gca, 'xlim', [0, 60]);
set(gca, 'ylim', [0, 180]);
set(gca, 'ytick', 0: 20: 180);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('Distance  (km)', 'Fontsize', 15);
text(2, 170, '(b)', 'fontsize', 18, 'color', 'k');
set(gca, 'position', [0.6 0.3 0.35 0.35]);

set(gcf, 'unit', 'centimeters', 'position', [5, 5, 20, 20], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/neic1d_hf.ps');
print('-depsc2','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/neic1d_hf.eps');
print('-dpdf','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/neic1d_hf.pdf');
% filename=['G:\BackProjection\mexico\AL-neic\abc.jpg' ];
% print('-djpeg', filename);
