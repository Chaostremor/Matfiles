% to mexico 3d high-freq. result, for new South America stations
% 
% Author: 
%     C. Song, 2018.11.15

%% initial setting
clear;
close all;
% % initial parameter, change when needed
% %%%%%%%%%%%%%%%%%%%
% lon0 = -93.9067;           % lat0, lon0 denote epicenter
% lat0 = 15.0356;
% data = load('/Users/chaos/PreviousResearch/BackProjection/mexico/NewSA_trans/mexico_nsa_neic453dstations10s0.5HzTo2Hz/HFdots_tc');          % HFdots dir
% %%%%%%%%%%%%%%%%%%%
% 
% %% determine the end of stacking time.
% tstack = data(:, end);
% power = data(:, 4);                    % normalized power, 0-1
% coeff =data(:, 5);
% alpha=power.*coeff;
% figure
% plot(tstack, alpha, 'linestyle', '-', 'color', 'k', 'linewidth', 2); hold on
% line([0 80], [0.1 0.1], 'linewidth', 0.5, 'color', [180/255 180/255 180/255], 'linestyle', '--'); hold on
% 
% %%
% ind=1:47;         % ending index, determined from alpha figure
% trup = data(ind, 1);                                 % time, 1st col.
% lat = data(ind, 2);                          % latitude, 2nd col.
% lon = data(ind, 3);                         % longitude, 2nd col.
% power = data(ind, 4);                    % normalized power, 0-1
% coeff =data(ind, 5);                       % cross-correlation coefficients, 0-1
% fdata = [trup lat lon power coeff];
% save('/Users/chaos/PreviousResearch/BackProjection/mexico/NewSA_trans/mexico_nsa_neic453dstations10s0.5HzTo2Hz/HFdots_final_nsa3d_hf', 'fdata','-ascii');
% save('/Users/chaos/PreviousResearch/BackProjection/mexico/NewSA_trans/mexico_nsa_neic453dstations10s0.5HzTo2Hz/nsa3dhf.mat');

load('/Users/chaos/PreviousResearch/BackProjection/mexico/NewSA_trans/nchile1.mat');
nm = ret.nm;
lat = ret.lat';
lon = ret.lon';
fid=fopen('/Users/chaos/PreviousResearch/BackProjection/mexico/NewSA_trans/main_nsa','w');
for i=1:length(lat)
% d(i, :) = strcat({nm(i,:)}, {num2str(lat(i))}, strcat{num2str(lon(i))});
d = [nm(i,:), num2str(lon(i)), '    ', num2str(lat(i))];
fprintf(fid,'%s \n',d);
end
% d