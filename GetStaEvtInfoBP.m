% get used stations an event info from alignteleBP result
% Author: C. Song, 2017.11.20

close all;

clear;
load('G:\BackProjection\mexico\SA_trans\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico\SA_trans\main_sa', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);

clear;
load('G:\BackProjection\mexico2\af1\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico2\af1\af1', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);

clear;
load('G:\BackProjection\mexico2\af2\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico2\af2\af2', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);

clear;
load('G:\BackProjection\mexico2\af3\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico2\af3\af3', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);

clear;
load('G:\BackProjection\mexico2\af4\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico2\af4\af4', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);

clear;
load('G:\BackProjection\mexico2\af5\nchile1.mat');
stnm=ret.nm;
lat=ret.lat;
lon=ret.lon;
fid=fopen('G:\BackProjection\mexico2\af5\af5', 'w');
for i = 1: length(lat)
    fprintf(fid, '%s %f     %f \n', stnm(i, :), lon(i), lat(i));
end
fclose(fid);