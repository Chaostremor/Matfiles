% This code is used to screen out stations from a set of stations
% station distribution may not be uniform, some part dense, while some
% sparse
% use distance and azimuth to choose one station within one small region
%
% Author: C. Song, 2017.9.17


clear; clc; 

% %%  way 1: select one ref. station, use dist and azi of others against it to restain
% NOT APPLICABLE: COST TOO MUCH MEMORY

% %1. get location of each station
% datadir = 'G:\BackProjection\mexico\EUO_trans\data\' ;    % 数据所在目录
% fid1 = fopen(strcat(datadir,'filelist')) ;      % strcat用于字符串连接
% file = textscan(fid1, '%s \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
% %
% fclose(fid1) ;
% fname = char(file{1});     % 台站名是第一列
% [s1, ~] = size(fname);
% %
% lon = zeros(s1,1);       % dist是头段中的震中距
% lat = zeros(s1,1);
% for ii = 1: s1
% % ii=1;
%     filename = strcat(datadir, fname(ii, :)) ;  % 获取文件路径
%     [~, ~, SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
%     lon(ii) = SAChdr.station.stlo ;   % 注意头段变量的存储方式
%     lat(ii) = SAChdr.station.stla;
% end
% % select a ref. station, the center of densest part
% lon0 = 8.151;
% lat0 = 46.76;
% [dist, azi] = distance (lat0, lon0, lat, lon);
% sdist = sort(dist);
% bd = 0.5;
% dd = 0.5;
% adist = bd: dd: max(dist);
% da = 90;
% s2 = length(adist);
% ind = [];
% % save('ssin.mat');
% for i = 1: s2
%     i
% %         i = 30;
%     %     indtemp = [];    
%     %     for j = 1: s1
%     %         if (dist(j) >= adist(i)-0.5*dd) && (dist(j)< adist(i)+0.5*dd)
%     %             indtemp = [indtemp; j];
%     %         end
%     %     end
%     temp1 = find((dist >= adist(i)-0.5*dd) & (dist< adist(i)+0.5*dd));
%     if length(temp1) > 1
%         azitemp = azi(temp1);
%         da = da*adist(1)/adist(i);
%         aazi = 0: da: 360-da;
%         s4 = length(aazi);
%         for j = 1: s4
% %             j
%             % j=4;
%             temp2 = find((azitemp >= aazi(j)-0.5*da) & (azitemp< aazi(j)+0.5*da));
%             if length(temp2) > 1
%                 [~, index] = min(abs(azitemp(temp2) - aazi(j)));
%                 ind = [ind; temp1(temp2(index))];
%             else
%                 ind = [ind; temp1(temp2)];
%             end
%             
%             %         [~, index] = min(abs(azitemp - aazi(j)));
%             %         ind = [ind; indtemp(index)];
%         end
%     else
%         ind = [ind; temp1];
%     end
% end
% 
% ind = unique(ind);
% fuse = fname(ind, :);
% latuse = lat(ind);
% lonuse = lon(ind);
% 
% figure
% plot(lon, lat, 'b.', 'MarkerSize', 20);
% axis equal
% 
% figure
% plot(lonuse, latuse,'b.', 'MarkerSize', 20);
% axis equal

%% way2: calculate dist against each station to others, and delete those dist smaller than threshold
%1. get location of each station
datadir = 'G:\BackProjection\mexico\EUOA_trans_disp\data\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'filelist1')) ;      % strcat用于字符串连接
file = textscan(fid1, '%s \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
fname = char(file{1});     % 台站名是第一列
[s1, ~] = size(fname);
%
lon = zeros(s1,1);       % dist是头段中的震中距
lat = zeros(s1,1);
for ii = 1: s1
% ii=1;
    filename = strcat(datadir, fname(ii, :)) ;  % 获取文件路径
    [~, ~, SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    lon(ii) = SAChdr.station.stlo ;   % 注意头段变量的存储方式
    lat(ii) = SAChdr.station.stla;
end
% use any station as ref.
mind = 0.5;
index = 1: s1;
index = index';
ind = [];
for i = 1: s1
    i    
    [dist, ~] = distance(lat(index(1)), lon(index(1)), lat(index(2: end)), lon(index(2: end)));
    temp = find(dist >= mind);
    if isempty(temp)
        ind = [ind; index(1)];
        break;
    elseif length(temp) == 1
        ind = [ind; index(1)];
        ind = [ind; index(temp+1)];
        break;
    else
        ind = [ind; index(1)];
        index = index(temp+1);
        continue;
    end    
end

fuse = fname(ind, :);
latuse = lat(ind);
lonuse = lon(ind);

figure
plot(lon, lat, 'b.', 'MarkerSize', 20);
axis equal

figure
plot(lonuse, latuse,'b.', 'MarkerSize', 20);
axis equal

fid2 = fopen(strcat(datadir,'filelist'), 'wt');
s2 = size(fuse);
for i = 1: s2
    fprintf(fid2, '%s \n', fuse(i, :));
end
fclose(fid2);
    

































