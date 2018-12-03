% check the cross-correlation procedure in cap
%
% Author: 
%     C. Song, 2017.7.20
% 

clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);

datadir = 'G:\Alxa\test2_13\alxa23_13_' ;    % 数据所在目录
t11 = zeros(175, sa);
t12 = zeros(175, sa);
data11 = zeros(175, sa);
data12 = zeros(175, sa);
t21 = zeros(175, sa);
t22 = zeros(175, sa);
data21 = zeros(175, sa);
data22 = zeros(175, sa);
npts = zeros(sa ,1);
for i =1:sa
% i=sa;
filename = strcat(datadir,strcat(stnm(i, :), '.8'));   % Pnl z, data
[t, data, ~] = fget_sac(filename) ;
npts(i) = length(data);
t11(1: npts(i), i) = t; 
data11(1: npts(i), i) = data;
filename = strcat(datadir,strcat(stnm(i, :), '.9'));   % Pnl z, syn
[t, data, ~] = fget_sac(filename) ;
t12(1: npts(i), i) = t; 
data12(1: npts(i), i) = data;

filename = strcat(datadir,strcat(stnm(i, :), '.6'));   % Pnl r, data
[t, data, ~] = fget_sac(filename) ;
npts(i) = length(data);
t21(1: npts(i), i) = t; 
data21(1: npts(i), i) = data;
filename = strcat(datadir,strcat(stnm(i, :), '.7'));   % Pnl r, syn
[t, data, ~] = fget_sac(filename) ;
t22(1: npts(i), i) = t; 
data22(1: npts(i), i) = data;
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

delta = 0.2;
% snlag = zeros(2*175+1, nuse);       % lag points 
% tsncorr = zeros(2*175+1, nuse);    % after correlation, length is 2*nsnlen+1
% dsncorr = zeros(2*175+1, nuse);

% deal with Pnl z 
t13 = zeros(175, sa);
for i = 1:sa
    t13(1: npts(i), i) = t12(1: npts(i), i)+tpshift(i);
end

t15 = zeros(175, sa);
for i = 1:sa
    t15(1: npts(i), i) = t12(1: npts(i), i)+totalshift(i);
end

i=24;
[coef1, lag1] = xcorr(data11(1: npts(i), i), data12(1: npts(i), i), 'coeff');
maxcoef1 = max(coef1);
index1 = find(coef1==maxcoef1);
nlag1 = lag1(index1);
tlag1 = nlag1*delta;
t14(1: npts(i), i) = t12(1: npts(i), i) + tlag1;
t16(1: npts(i), i) = t12(1: npts(i), i) + tlag1 + tpshift(i);

figure
plot(t11(1: npts(i), i), data11(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'g'); hold on             % syn data
% plot(t3(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on              % syn data shifted by tpshift
plot(t14(1: npts(i), i), data12(1: npts(i), i), 'LineWidth',2, 'color', 'b'); hold on             % syn data shifted by crscrl lag time 
plot(t15(1: npts(i), i), data12(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             %  syn data shifted by CAP (cap's crscrl + tpshift)
plot(t16(1: npts(i), i), data12(1: npts(i), i), 'LineWidth',2, 'color', 'm'); hold on            % syn data shifted by crscrl lag time + tpshift
%%% it turns out that "r" and "m" is the same process

% deal with Pnl r
t23 = zeros(175, sa);
for i = 1:sa
    t23(1: npts(i), i) = t22(1: npts(i), i)+tpshift(i);
end

t25 = zeros(175, sa);
for i = 1:sa
    t25(1: npts(i), i) = t22(1: npts(i), i)+totalshift(i);
end

i=24;
[coef2, lag2] = xcorr(data21(1: npts(i), i), data22(1: npts(i), i), 'coeff');
maxcoef2 = max(coef2);
index2 = find(coef2==maxcoef2);
nlag2 = lag2(index2);
tlag2 = nlag2*delta;
t24(1: npts(i), i) = t22(1: npts(i), i) + tlag2;
t26(1: npts(i), i) = t22(1: npts(i), i) + tlag2 + tpshift(i);

figure
plot(t21(1: npts(i), i), data21(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'g'); hold on             % syn data
% plot(t3(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on              % syn data shifted by tpshift
plot(t24(1: npts(i), i), data22(1: npts(i), i), 'LineWidth',2, 'color', 'b'); hold on             % syn data shifted by crscrl lag time 
plot(t25(1: npts(i), i), data22(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             %  syn data shifted by CAP (cap's crscrl + tpshift)
plot(t26(1: npts(i), i), data22(1: npts(i), i), 'LineWidth',2, 'color', 'm'); hold on            % syn data shifted by crscrl lag time + tpshift