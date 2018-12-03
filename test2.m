% fit model check
%
% USAGE:
%     check if the result model from travel_time_curve_fitting is reasonable,    
%     plot the travel time curve on the data based on the new model
%
%  Author:   C. Song, 2017.5.18
%       

%% 1. read data
datadir = 'G:\Alxa\real\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:), '.z')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t1 = zeros(sc,sa);          % t是时间，data是数据
data1 = zeros(sc,sa);   % data1 存放Z分量数据
data2 = zeros(sc, sa);  % data2 存放T分量数据
dist = zeros(sa,1);       % dist是头段中的震中距
pickpg = zeros(sa,1);        % pg是直达P波到时
picksg = zeros(sa,1);         % sg是直达S波到时
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:), '.z')) ;  % 获取z分量数据路径
    [t1(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data1是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    az(i) = SAChdr.evsta.az;
    stla(i) = SAChdr.station.stla;
    stlo(i) = SAChdr.station.stlo;    
    pickpg(i) = SAChdr.times.t2;
    picksg(i) = SAChdr.times.t3;
    filename = strcat(datadir,strcat(stnm(i,:), '.t')) ;  % 获取t分量数据路径
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end
evla = SAChdr.event.evla;
evlo = SAChdr.event.evlo;
% fid2 = fopen('G:\Alxa\nodecimate\3test400\timemark3.dat');
% time = textscan(fid2, '%s %f %f %f \n');
% fclose(fid2) ;
% cappg = time{4};
% shift = time{3};
fid2 = fopen('G:\Alxa\nodecimate\3test400\timeshiftfromcap');
time = textscan(fid2, '%s %f\n');
fclose(fid2) ;
shift = time{2};

datadir = 'G:\Alxa\syn\' ;    % 数据所在目录
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.z'));  % 获取z分量数据路径
    [t3(:, i), data3(:, i), ~] = fget_sac(filename); 
end

vs = [3.34; 3.61; 4.11; 4.48];
vp = [5.2; 6.2; 6.4; 8.12];
tk= [7.7; 30.3; 16.2; 0.0];
dep = 18;

%% 3. calculate P wave phase travel time curve
distsamp = 0: 0.1: 500;

% for direct P, Pg, model 
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
pgint = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve

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
sgint = interp1(xsg, tsg, dist, 'spline');

figure
for i = 1:sa
%     plot(t1(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
%     plot(t3(:,i), data3(:,i)./(max(abs(data3(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'r', 'LineWidth',1);hold on;
    plot(t1(:,i), data1(:,i) .*1000.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
    plot(t3(:,i), data3(:,i) .*1000.0 + double(dist(i)), 'linestyle', '-', 'color', 'r', 'LineWidth',1);hold on;
%     plot(cappg(i), dist(i), 'g.', 'MarkerSize', 20); hold on;           % pg picks
end
plot(tpg, xpg, 'linestyle', '--', 'color', 'b', 'LineWidth', 2); hold on;                      % pg model
% plot(pickpg, dist, 'g.', 'MarkerSize', 8); hold on;                                                    % pg picks                                         
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 20: 420);

figure
plot(dist, shift, 'k*', 'MarkerSize', 10);

% figure
% plot(az, shift, 'k*', 'MarkerSize', 10);
% 
% figure
% plot(evlo, evla,'ko', 'MarkerSize', 10); hold on
% for i=1: sa
%     if shift(i) < 0
%         plot(stlo(i), stla(i),'r*', 'MarkerSize', 10 ); hold on
%     else
%         plot(stlo(i), stla(i),'b*', 'MarkerSize', 10 ); hold on
%     end
% end