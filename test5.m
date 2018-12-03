% test the relationship between real and synthetic data from cap
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

fid2 = fopen('G:\Alxa\test2_13\timemark5.dat') ;      % strcat用于字符串连接
time = textscan(fid2, '%s %f %f %f %f %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid2) ;
synpg = time{2};
synsg = time{3};
tpshift = time{4};     % 台站名是第一列
corshift = time{5};
totalshift = time{6};

% fid3 = fopen('G:\Alxa\nodecimate\3test400\timemark4.dat') ;      % strcat用于字符串连接
% time = textscan(fid3, '%s %f %f %f  %f %f\n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
% %
% fclose(fid3) ;
% synpg = time{2};
% synsg = time{3};

datapg = synpg+totalshift;
datasg = synsg+totalshift;

t3 = zeros(175, sa);
for i = 1:sa
    t3(1: npts(i), i) = t2(1: npts(i), i)+totalshift(i);
end


% Pnl (r)
for i =1:sa
    figure (2)
% i = 2;
i
% ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
    plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on
    plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on
    plot(synpg(i), 0, 'b.','MarkerSize', 40);
%     plot(t3(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'g'); hold on
%     plot(datapg(i), 0, 'b.','MarkerSize', 40);
    pause (5);
    close (2)
end
