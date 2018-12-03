% plot mexico 2017 high-freq nomalized power, version 3
% with 1d
%
% Author: C. Song, 2018.11.22

clear;
close all;

% load result file
a = load('/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_al_neic210stations10s0.5HzTo2Hz/HFdots_tc');
t1d = a(:, 6);
p1d = a(:, 4);
c1d = a(:, 5);
nt = 1: 0.001: t1d(end);
alp1d = p1d.*c1d;
nalp1d = interp1(t1d, alp1d, nt, 'spline');

figure
% l1 = plot(t1d, alp1d, 'linestyle', '-', 'color', [130/255 130/255 130/255], 'linewidth', 2); hold on
% l2 = plot(t3d, alp3d, 'linestyle', '-', 'color', 'k', 'linewidth', 2); hold on
% l1 = plot(nt, nalp1d, 'linestyle', '-', 'color', [130/255 130/255 130/255], 'linewidth', 2); hold on
l1 = plot(nt, nalp1d, 'linestyle', '-', 'color', 'k', 'linewidth', 2); hold on
line([47.5 47.5], [0 0.8], 'linewidth', 1, 'color', [180/255 180/255 180/255], 'linestyle', '--'); hold on
line([0 80], [0.1 0.1], 'linewidth', 1, 'color', [180/255 180/255 180/255], 'linestyle', '--'); hold on
arrow([32, 0.2], [32, 0.32], 10, 'Color', [180/255 180/255 180/255], 'Width', 1.5, 'BaseAngle', 60); hold on
set(gca, 'fontsize', 12);
set(gca, 'XTick', 0: 10: 80);
set(gca, 'xlim', [0, 80]);
set(gca, 'ylim', [0, 0.8]);
set(gca, 'ytick', 0: 0.1: 0.8);
xlabel('Time  (s) ', 'Fontsize', 15);
ylabel('p*c', 'Fontsize', 15);
% legend([l1, l2],'1-D','3-D', 'location','northeast', 'fontsize', 20);
% hl=legend('1-D','3-D');
hl=legend('1-D');
set(hl, 'location','northeast', 'fontsize', 13);
text(15, 0.15, 'Stage 1', 'fontsize', 12, 'color', 'k');
text(35, 0.15, 'Stage 2', 'fontsize', 12, 'color', 'k');


set(gcf, 'unit', 'centimeters', 'position', [5, 5, 16, 10], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_power1d.ps');
print('-depsc2','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_power1d.eps');
print('-dpdf','-r600','/Users/chaos/PreviousResearch/BackProjection/mexico/AL-neic/mexico_power1d.pdf');