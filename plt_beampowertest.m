% plot beamforming power get by runteleBP, referring to movieteleBP
% Author: C. Song, 2017.12.3
%

clear;
% close all;


cd('G:\BackProjection\mexico\deptest2\mexico_neic_depth2103dstations10s0.5HzTo2Hz');
figure
row = 2;
col = 3;
step = 10;
tstart= 10;
tmark=45;               % specific time window
for i = 1: row
    for j = 1: col
        num = (i-1)*col+j;              % subplot number
        kt = tstart+(num-1)*step;   % source depth
        subplot(row, col, num);
        [pow]=plt_BeamPowtest(kt, num, tmark);
        set(gca, 'XTick', -95: 1: -92);
        set(gca, 'ytick', 14: 1: 17);
        set(gca, 'Fontsize', 12);
        xlabel('Longitude (бу)', 'Fontsize', 15);
        ylabel('Latitude (бу)', 'Fontsize', 15);
    end
end
set(gcf, 'unit', 'centimeters', 'position', [12, 5, 25, 15], 'PaperPositionMode', 'auto');
% print('-dpsc','-r600','G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\beamsnapshot.ps');
% print('-dpdf','-r600','G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\beamsnapshot.pdf');

print('-dtiff','-r300','G:\BackProjection\mexico\deptest2\mexico_neic_depth2103dstations10s0.5HzTo2Hz\depthsnapshot.tif');
print('-depsc2','-r600','G:\BackProjection\mexico\deptest2\mexico_neic_depth2103dstations10s0.5HzTo2Hz\depthsnapshot.eps');
print('-dpdf','-r600','G:\BackProjection\mexico\deptest2\mexico_neic_depth2103dstations10s0.5HzTo2Hz\depthsnapshot.pdf');