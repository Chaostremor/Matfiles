% plot real and synthetic data
% 
% USAGE:
%     plot real and synthetic data, 
% 
% Author: C. Song, 2017.5.31

%% 1. read data
realdir = 'G:\Alxa\realsynalxa39filter0.05to0.2\' ;    % 数据所在目录
% fid1 = fopen(strcat(realdir,'nweight.dat')) ;      % strcat用于字符串连接
fid1 = fopen(strcat(realdir,'nweight_selected.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, ~] = size(stnm);
% test to get size of data
ii=1;
filename = strcat(realdir,strcat(stnm(ii,:),'.t')) ;  % 获取文件路径
[ttest, ~, SAChdrtest] = fget_sac(filename) ;
dt1 = SAChdrtest.times.delta ;
npts1 = SAChdrtest.data.trcLen ; 
%
[sc1, ~] = size(ttest);
t1 = zeros(sc1, sa);          % t是时间，data是数据
data1 = zeros(sc1, sa);
dist = zeros(sa, 1);       % dist是头段中的震中距
for ii = 1: sa
    filename = strcat(realdir,strcat(stnm(ii,:),'.t')) ;  % 获取文件路径
    [t1(:,ii), data1(:,ii), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(ii) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end

syndir = 'G:\Alxa\realsynalxa39filter0.05to0.2\' ;    % 数据所在目录
% test to get size of data
ii=1;
filename = strcat(syndir,strcat(stnm(ii,:),'_syn.t')) ;  % 获取文件路径
[ttest, ~, SAChdrtest] = fget_sac(filename) ;
dt2 = SAChdrtest.times.delta ;
npts2 = SAChdrtest.data.trcLen ; 
%
[sc2, ~] = size(ttest);
t2 = zeros(sc2, sa);          % t是时间，data是数据
data2 = zeros(sc2, sa);
for ii = 1: sa
    filename = strcat(syndir,strcat(stnm(ii,:),'_syn.t')) ;  % 获取文件路径
    [t2(:,ii), data2(:,ii), ~] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
end

%% 2. plot data, real and syn.
%%% no generalition, multify 1000 to amplify, can see attenuation in real
%%% and syn data
% figure
% for i = 1:sa
% % i=1;
% %     data2(:, i) = data2(:, i).*3.*(dist(i)/100)^1; 
%     plot(t1(:,i), data1(:,i).*1000 + double(dist(i)), 'k-', 'LineWidth',1); hold on;
%     plot(t2(:,i), data2(:,i).*1000 + double(dist(i)), 'r-', 'LineWidth',1); hold on;
% end
% set(gca, 'XLIM', [0, 130]);
% set(gca, 'XTICK', 0: 10: 130);

%% 3. check the correctness of travel time in synthetic data
% model para.
vs = [2.09; 3.56; 4.12; 4.48];
vp = [3.70; 5.94; 6.30; 8.12];
tk= [6.1; 25.2; 14.9; 0.0];
dep = 13;

% syn Pg time 
nsample1 = 10000;
rayp1 = linspace(0, 1.0/vp(2)-0.0000005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for ii = 1: nsample1
    etap1 = sqrt((1.0/vp(1))^2-rayp1(ii)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(ii)^2);
    xpg(ii) = rayp1(ii) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
    tpg(ii) = rayp1(ii)*xpg(ii) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
synpg = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve 

% for Pn, model
nsample2 = 1000;
rayp2 =  1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
dxpn = linspace(0, 400, nsample2)';
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
synpn = interp1(xpn, tpn, dist, 'spline'); 

% for pPn, model
xppnmin = rayp2 * ( 3*tk(1)/etap1 + (dep-tk(1)+2*tk(2))/etap2 + 2*tk(3)/etap3 );
tppnmin = rayp2 * xppnmin +3*tk(1)*etap1 + (dep-tk(1)+2*tk(2))*etap2 + 2*tk(3)*etap3;
dxppn = linspace(0, 400, nsample2)';
xppn = xppnmin +dxppn;
tppn = tppnmin +rayp2 .*dxppn;
synppn = interp1(xppn, tppn, dist, 'spline');

% PmP, moho reflection
rayp3 = linspace(0, 1.0/vp(3)-0.0000005, nsample1)';
xpmp = zeros(nsample1,1);
tpmp = zeros(nsample1,1);
for ii = 1: nsample1
    etap1 = sqrt((1.0/vp(1))^2-rayp3(ii)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp3(ii)^2);
    etap3 = sqrt((1.0/vp(3))^2-rayp3(ii)^2);
    xpmp(ii) = rayp3(ii) * ( tk(1)/etap1 + tk(2)/etap2+ 2*tk(3)/etap3 + (tk(2)+tk(1)-dep)/etap2 );
    tpmp(ii) = rayp3(ii)*xpmp(ii) + tk(1)*etap1+ tk(2)*etap2+ 2*tk(3)*etap3 + (tk(2)+tk(1)-dep)*etap2;
end
synpmp = interp1(xpmp, tpmp, dist, 'spline');       % model determined pg curve 

% for direct S, Sg
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for ii = 1: nsample1
    etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
    xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
    tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
synsg = interp1(xsg, tsg, dist, 'spline');

% for Sn
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
dxsn = linspace(0, 400, nsample2)';
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;           
synsn = interp1(xsn, tsn, dist, 'spline');       % model determined sn curve

syndir = 'G:\Alxa\synalxa23nofilter\' ;    % 数据所在目录
% test to get size of data
ii=1;
filename = strcat(syndir,strcat(stnm(ii,:),'_syn.z')) ;  % 获取文件路径
[ttest, ~, SAChdrtest] = fget_sac(filename) ;
dt3 = SAChdrtest.times.delta ;
npts3 = SAChdrtest.data.trcLen ; 
%
[sc3, ~] = size(ttest);
t3 = zeros(sc3, sa);          % t是时间，data是数据
data3 = zeros(sc3, sa);
for ii = 1: sa
    filename = strcat(syndir,strcat(stnm(ii,:),'_syn.z')) ;  % 获取文件路径
    [t3(:,ii), data3(:,ii), ~] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
end

%%%  with generalition
figure
for i = 1:sa
% i=1;
%     plot(t2(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'r-', 'LineWidth',1); hold on;
    plot(t3(:,i), data3(:,i)./(max(abs(data3(:,i)))) .*15.0 + double(dist(i)), 'r-', 'LineWidth',1); hold on;
    plot(synpg, dist, 'k.', 'MarkerSize', 8); hold on;
    plot(synpn, dist, 'b.', 'MarkerSize', 8); hold on;
    plot(synppn, dist, 'y.', 'MarkerSize', 8); hold on;
%     plot(synpmp, dist, 'y.', 'MarkerSize', 8); hold on;
%     plot(synsg, dist, 'k.', 'MarkerSize', 8); hold on;
%     plot(synsn, dist, 'b.', 'MarkerSize', 8); hold on;
end
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 10: 130);




























