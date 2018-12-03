% simulate a model for 3D mexico BP
%
% C. Song, 2018.3.20
clear; close all;
load('G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\al3dhf.mat');
t1 = mean(tuse(1:7));
lat1= mean(latuse(1:7));
lon1 = mean(lonuse(1:7));
pow1 = mean(powuse(1:7));
t2 = mean(tuse(8:15));
lat2 = mean(latuse(8:15));
lon2 = mean(lonuse(8:15));
pow2 = mean(powuse(8:15));
t3 = mean(tuse(16:23));
lon3 = mean(lonuse(16:23));
lat3 = mean(latuse(16:23));
pow3 = mean(powuse(16:23));
t4 = mean(tuse(24: 31));
lon4 = mean(lonuse(24: 31));
lat4 = mean(latuse(24: 31));
pow4 = mean(powuse(24: 31));
t5 = mean(tuse(32: 35));
lon5 = mean(lonuse(32: 35));
lat5 = mean(latuse(32: 35));
pow5 = mean(powuse(32: 35));
t6 = mean(tuse(36: 41));
lon6 = mean(lonuse(36: 41));
lat6 = mean(latuse(36: 41));
pow6 = mean(powuse(36: 41));
barrier = [t1 lon1 lat1 pow1; t2 lon2 lat2 pow2; t3 lon3 lat3 pow3; t4 lon4 lat4 pow4; t5 lon5 lat5 pow5; t6 lon6 lat6 pow6];
save('G:\BackProjection\mexico\AL\mexico_al2103dstations10s0.5HzTo2Hz\barrier', 'barrier','-ascii');