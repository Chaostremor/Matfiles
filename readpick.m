% read time mark of data, to fit the picks to a curve
%

clear;
fid1 = fopen('L:\nodecimate\400nofManPick\timemarknofnod.dat') ;      % strcat用于字符串连接
time = textscan(fid1, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid1) ;
stnm = char(time{1});
dist = time{2};
realpg = time{3};
realsg = time{4};

% Pg picks fitting curve 
distsamp = 0: 0.1: 500;
a1 = 6.104;
b1 = 36.1;
% a1 = 6.443;
% b1 = 39.43;
pgfit = a1.*sqrt((distsamp./b1).^2+1);
raypfit = a1/b1.*distsamp./sqrt(distsamp.^2+b1^2);

a2 = 10.53;
b2 = 37.2;
sgfit = a2.*sqrt((distsamp./b2).^2+1);
raysfit = a2/b2.*distsamp./sqrt(distsamp.^2+b2^2);

% %% 1. read data, lowpass 0.5 hz
% datadir = 'J:\台式机G盘\Alxa\nodecimate\3test400\' ;    % 数据所在目录
% % datadir = 'G:\Alxa\nodecimate\3test400_mis_alxa23\' ;    % 数据所在目录
% fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
% weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
% %
% fclose(fid1) ;
% stnm = char(weight{1});     % 台站名是第一列
% [sa, sb] = size(stnm);
% % test to get size of data
% i=1;
% filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
% [ttest, datatest, SAChdrtest] = fget_sac(filename) ;
% disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
% dt = SAChdrtest.times.delta ;
% npts = SAChdrtest.data.trcLen ; 
% %
% [sc, sd] = size(ttest);
% t = zeros(sc,sa);          % t是时间，data是数据
% data1 = zeros(sc,sa);   % data1 存放Z分量数据
% data2 = zeros(sc, sa);  % data2 存放T分量数据
% dist = zeros(sa,1);       % dist是头段中的震中距
% pg = zeros(sa,1);        % pg是直达P波到时
% sg = zeros(sa,1);         % sg是直达S波到时
% for i = 1:sa
%     filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取z分量数据路径
%     [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data1是数据，SAChdr是头段变量
%     dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
%     pg(i) = SAChdr.times.t2;
%     sg(i) = SAChdr.times.t3;
%     filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取t分量数据路径
%     [~, data2(:,i), ~] = fget_sac(filename) ; 
% end
% 
% figure
% for i = 1:sa
%     plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % pg fit
% set(gca, 'XLIM', [0, 100]);
% set(gca, 'XTICK', 0: 10: 100);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);

% figure
% for i = 1:sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % sg fit
% set(gca, 'XLIM', [0, 100]);
% set(gca, 'XTICK', 0: 10: 100);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);

%% 2. read data, bandpass 0.05-0.2 hz
datadir = 'J:\台式机G盘\Alxa\realsynalxa39filter0.05to0.2\' ;    % 数据所在目录
% datadir = 'G:\Alxa\nodecimate\3test400_mis_alxa23\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data1 = zeros(sc,sa);   % data1 存放Z分量数据
data2 = zeros(sc, sa);  % data2 存放T分量数据
dist = zeros(sa,1);       % dist是头段中的震中距
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取z分量数据路径
    [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data1是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取t分量数据路径
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end

figure
for i = 1:sa
    plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
plot(pgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 20: 420);

% figure
% for i = 1:sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % sg fit
% set(gca, 'XLIM', [0, 150]);
% set(gca, 'XTICK', 0: 10: 150);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);