% plot beamforming power get by runteleBP, referring to movieteleBP
% Author: C. Song, 2017.12.3
%

clear;
close all;

% % way 1 to save fig, imwrite
% h=figure;
% kt=44;
% plt_BeamPow(kt);
% colorbar;
% set(gca, 'XTick', -95: 1: -92);
% set(gca, 'ytick', 14: 1: 17);
% set(gca, 'Fontsize', 12);
% xlabel('Longitude ^o', 'Fontsize', 15);
% ylabel('Latitude ^o', 'Fontsize', 15);
% set(gcf, 'unit', 'centimeters', 'position', [15, 10, 12, 10], 'PaperPositionMode', 'auto');
% f = getframe(h);
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% imwrite(im, map, 'G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\aa.tiff', 'tiff', 'Resolution', 600);

% % way 2 to save fig, print
% figure;
% kt=45;
% [pow]=plt_BeamPow(kt);
% minp=min(pow(:));
% maxp=max(pow(:));
% colorbar;
% % c=colormap;
% caxis([minp, maxp]);
% set(gca, 'XTick', -95: 1: -92);
% set(gca, 'ytick', 14: 1: 17);
% set(gca, 'Fontsize', 12);
% xlabel('Longitude (бу)', 'Fontsize', 15);
% ylabel('Latitude (бу)', 'Fontsize', 15);
% set(gcf, 'unit', 'centimeters', 'position', [15, 10, 12, 10], 'PaperPositionMode', 'auto');
% % print('-dpsc','-r600','G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\beampower.ps');
% % print('-dpdf','-r600','G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\beampower.pdf');
% print('-dtiff','-r300','G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\beampower.tif');

% % way 3 to save fig, saveas, cannot specialize resolution
% h=figure;
% kt=44;
% plt_BeamPow(kt);
% colorbar;
% set(gca, 'XTick', -95: 1: -92);
% set(gca, 'ytick', 14: 1: 17);
% set(gca, 'Fontsize', 12);
% xlabel('Longitude ^o', 'Fontsize', 15);
% ylabel('Latitude ^o', 'Fontsize', 15);
% set(gcf, 'unit', 'centimeters', 'position', [15, 10, 12, 10], 'PaperPositionMode', 'auto');
% saveas(h, 'G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\cc.tif', 'tiff');
% saveas(h, 'G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\cc.pdf', 'pdf');
% export_fig G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\aaa.pdf -pdf -depsc -transparent -r100

cd('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz');
figure
row = 3;
col = 3;
step = 5;
tstart= 5;
for i = 1: row
    for j = 1: col
        num = (i-1)*col+j;
        kt = tstart+(num-1)*step;
        subplot(row, col, num);
        [pow]=plt_BeamPow(kt);
        set(gca, 'XTick', -95: 1: -92);
        set(gca, 'ytick', 14: 1: 17);
        set(gca, 'Fontsize', 12);
        xlabel('Longitude (бу)', 'Fontsize', 15);
        ylabel('Latitude (бу)', 'Fontsize', 15);
    end
end
set(gcf, 'unit', 'centimeters', 'position', [12, 5, 25, 23], 'PaperPositionMode', 'auto');
print('-depsc2','-r600','G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\beamsnapshot.eps');
print('-dpdf','-r600','G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\beamsnapshot.pdf');
print('-dtiff','-r300','G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\beamsnapshot.tif');