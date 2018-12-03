% test 6
% check relationship between the syn data cut by CAP and arrival time calc. by matlab
%

clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\test2_13\alxa23_13_' ;    % 数据所在目录
t1 = zeros(175, sa);
t2 = zeros(175, sa);
data1 = zeros(175, sa);
data2 = zeros(175, sa);
npts = zeros(sa ,1);
for i =1:sa
% i=sa;
filename = strcat(datadir,strcat(stnm(i, :), '.8'));   % Pnl z, data
[t, data, ~] = fget_sac(filename) ;
npts(i) = length(data);
t1(1: npts(i), i) = t; 
data1(1: npts(i), i) = data;
filename = strcat(datadir,strcat(stnm(i, :), '.9'));   % Pnl z, syn
[t, data, ~] = fget_sac(filename) ;
t2(1: npts(i), i) = t; 
data2(1: npts(i), i) = data;
end
dt = t2(2, 1) - t2(1, 1);

fid2 = fopen('G:\Alxa\test2_13\timemark5.dat') ;      % strcat用于字符串连接
time = textscan(fid2, '%s %f %f %f %f %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid2) ;
synpg = time{2};
synsg = time{3};

figure
for i = 1: 3: sa
% i=sa;
    plot(t1(1: npts(i), i), data1(1: npts(i), i)./max(data1(1: npts(i), i))*15+double(dist(i)), 'LineWidth',2, 'color', 'k'); hold on
    plot(t2(1: npts(i), i), data2(1: npts(i), i)./max(data2(1: npts(i), i))*15+double(dist(i)), 'LineWidth',2, 'color', 'r'); hold on             % syn data
    plot(synpg(i), double(dist(i)), 'b.','MarkerSize', 40);
    plot(synsg(i), double(dist(i)), 'g.','MarkerSize', 40);
end

t3 = zeros(400, sa);
t4 = zeros(400, sa);
data3 = zeros(400, sa);
data4 = zeros(400, sa);
npts = zeros(sa ,1);
for i =1:sa
% i=sa;
filename = strcat(datadir,strcat(stnm(i, :), '.0'));   % Pnl z, data
[t, data, ~] = fget_sac(filename) ;
npts(i) = length(data);
t3(1: npts(i), i) = t; 
data3(1: npts(i), i) = data;
filename = strcat(datadir,strcat(stnm(i, :), '.1'));   % Pnl z, syn
[t, data, ~] = fget_sac(filename) ;
t4(1: npts(i), i) = t; 
data4(1: npts(i), i) = data;
end
dt = t4(2, 1) - t4(1, 1);

fid2 = fopen('G:\Alxa\test2_13\timemark5.dat') ;      % strcat用于字符串连接
time = textscan(fid2, '%s %f %f %f %f %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid2) ;
synpg = time{2};
synsg = time{3};

figure
for i = 1: 3: sa
% i=sa;
    plot(t3(1: npts(i), i), data3(1: npts(i), i)./max(data3(1: npts(i), i))*15+double(dist(i)), 'LineWidth',2, 'color', 'k'); hold on
    plot(t4(1: npts(i), i), data4(1: npts(i), i)./max(data4(1: npts(i), i))*15+double(dist(i)), 'LineWidth',2, 'color', 'r'); hold on             % syn data
    plot(synsg(i), double(dist(i)), 'b.','MarkerSize', 40);
end