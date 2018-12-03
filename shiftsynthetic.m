clear ; clc ;
close all;
%% read data, real or synthetic

datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为weight{1}  
fclose(fid1) ;
stnm = char(weight{1});
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.r')) ;  % 获取文件路径
%filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
% load data
[sc, sd] = size(ttest);
treal = zeros(sc,sa);
real = zeros(sc,sa);
dist = zeros(sa,1);
azi = zeros(sa,1);
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.r')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [treal(:,i), real(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    azi(i) = SAChdr.evsta.az ;
end

datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录
%tsyn = zeros(sc,sa);
%syn = zeros(sc,sa);
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.r')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [tsyn(:,i), syn(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
end

fid1 = fopen('G:\Alxa\tshift.txt') ;      % strcat用于字符串连接
tshift = textscan(fid1, '%s %f %f %f %f %f %f %f') ;  % 得到的是一个cell array,调用语句为weight{1}  
fclose(fid1) ;
pzshift = tshift{4};
prshift = tshift{5};
for i = 1:sa
    tssyn(:,i) = tsyn(:,i)+prshift(i);
end
figure (1)
for i = 1:sa
    plot(treal(:,i), real(:,i)./(max(abs(real(:,i)))) .*15.0 + double(dist(i)), 'black');hold on;
    plot(tsyn(:,i), syn(:,i)./(max(abs(syn(:,i)))) .*15.0 + double(dist(i)), 'red');hold on;
    plot(tssyn(:,i), syn(:,i)./(max(abs(syn(:,i)))) .*15.0 + double(dist(i)), 'green');hold on;
end
set(gca,'YLIM', [round(min(dist))-10, round(max(dist))+10]);
set(gca,'YLIM', [-10, 500]);