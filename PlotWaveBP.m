% plot used stations' seismograms from alignteleBP
% Author: C. Song, 2017.11.20

close all;

%% mainshock,  AL
clear;
load('G:\BackProjection\mexico\AL\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
line([50 50], [-10 220], 'linewidth', 2.5, 'color', 'k', 'linestyle', '--'); hold on
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [0, 80]);
set(gca, 'xtick', 0: 10: 80);
set(gca, 'ylim', [-10, 220]);
set(gca, 'ytick', 0: 40: 200);
text(1, 215, '(a)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 18], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico\AL\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico\AL\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico\AL\waveform.pdf');
% figure
% imagesc()

%% mainshock, NA
% clear;
% load('G:\BackProjection\mexico\SA_trans\nchile1.mat');
% ret.plotcolor=[70/255 70/255 70/255];
% figure
% plotwave(ret);
% line([50 50], [-10 50], 'linewidth', 2.5, 'color', 'k', 'linestyle', '--'); hold on
% xlabel('Time (s)', 'Fontsize', 18);
% ylabel('Trace index', 'Fontsize', 18);
% set(gca,'FontSize',15)
% set(gca, 'xlim', [0, 80]);
% set(gca, 'xtick', 0: 10: 80);
% set(gca, 'ylim', [-10, 50]);
% set(gca, 'ytick', 0: 20: 50);
% set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 6], 'PaperPositionMode', 'auto');
% print('-dpsc','-r600','G:\BackProjection\mexico\SA_trans\waveform.ps');
% print('-dpdf','-r600','G:\BackProjection\mexico\SA_trans\waveform.pdf');



%% aftershocks, af1
clear;
load('G:\BackProjection\mexico2\af1\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [-5, 70]);
set(gca, 'xtick', 0: 10: 70);
set(gca, 'ylim', [-10, 80]);
set(gca, 'ytick', 0: 20: 80);
text(-4, 76, '(a)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 10], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico2\af1\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico2\af1\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico2\af1\waveform.pdf');

%% aftershocks, af2
clear;
load('G:\BackProjection\mexico2\af2\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [-5, 70]);
set(gca, 'xtick', 0: 10: 70);
set(gca, 'ylim', [-10, 165]);
set(gca, 'ytick', 0: 20: 165);
text(-4, 161, '(b)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 18], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico2\af2\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico2\af2\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico2\af2\waveform.pdf');

%% aftershocks, af3
clear;
load('G:\BackProjection\mexico2\af3\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [-5, 70]);
set(gca, 'xtick', 0: 10: 70);
set(gca, 'ylim', [-10, 135]);
set(gca, 'ytick', 0: 20: 135);
text(-4, 131, '(c)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 15], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico2\af3\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico2\af3\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico2\af3\waveform.pdf');

%% aftershocks, af4
clear;
load('G:\BackProjection\mexico2\af4\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [-5, 70]);
set(gca, 'xtick', 0: 10: 70);
set(gca, 'ylim', [-10, 65]);
set(gca, 'ytick', 0: 20: 65);
text(-4, 61, '(d)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 8.5], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico2\af4\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico2\af4\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico2\af4\waveform.pdf');

%% aftershocks, af5
clear;
load('G:\BackProjection\mexico2\af5\nchile1.mat');
ret.plotcolor=[70/255 70/255 70/255];
figure
plotwave(ret);
xlabel('Time (s)', 'Fontsize', 18);
ylabel('Trace index', 'Fontsize', 18);
set(gca,'FontSize',15)
set(gca, 'xlim', [-5, 70]);
set(gca, 'xtick', 0: 10: 70);
set(gca, 'ylim', [-10, 60]);
set(gca, 'ytick', 0: 20: 60);
text(-4, 56, '(e)', 'fontsize', 18, 'color', 'k');
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 24, 8], 'PaperPositionMode', 'auto');
print('-dpsc2','-r600','G:\BackProjection\mexico2\af5\waveform.ps');
print('-depsc2','-r600','G:\BackProjection\mexico2\af5\waveform.eps');
print('-dpdf','-r600','G:\BackProjection\mexico2\af5\waveform.pdf');

