% fit model check
%
% USAGE:
%     check if the result model from travel_time_curve_fitting is reasonable,    
%     plot the travel time curve on the data based on the new model
%
%  Author:   C. Song, 2017.5.18
%       

%% 1. read data
datadir = 'G:\Alxa\nodecimate\3test400\' ;    % 数据所在目录
% datadir = 'G:\Alxa\nodecimate\3test400_reloc\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
% fid1 = fopen(strcat(datadir,'reweight.dat')) ;      % strcat用于字符串连接
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
pg = zeros(sa,1);        % pg是直达P波到时
sg = zeros(sa,1);         % sg是直达S波到时
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取z分量数据路径
    [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data1是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    pg(i) = SAChdr.times.t2;
    sg(i) = SAChdr.times.t3;
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取t分量数据路径
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end
fid2 = fopen('G:\Alxa\nodecimate\3test400\timemark2.dat') ;      % strcat用于字符串连接
time = textscan(fid2, '%s %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid2) ;
pgo = time{2};     % 台站名是第一列
sgo = time{3};

% fid3 = fopen('G:\Alxa\nodecimate\3test400\timemark5.dat') ;      % strcat用于字符串连接
% time = textscan(fid3, '%s %f %f %f %f %f %f %f\n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
% %
% fclose(fid3) ;
% tpshift = time{4};
% corshift = time{5};
% totalshift = time{6};     % 台站名是第一列

fid4 = fopen('G:\Alxa\nodecimate\3test400_reloc\timeshiftfromcap32_18') ;      % strcat用于字符串连接
time = textscan(fid4, '%s %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid4) ;
tpshift2 = time{2};
corshift2 = time{3};     % 台站名是第一列
totalshift2 = tpshift2 + corshift2;

% distuse = [];
% shiftuse = [];
% pguse = [];
% sguse = [];
% stnmuse = [];
% for i = 1:sa
%     if totalshift(i) > -5 && totalshift(i) < 1
%         distuse =  [distuse; dist(i)];
%         shiftuse = [shiftuse; totalshift(i)];
%         pguse = [pguse; pg(i)];
%         sguse = [sguse; sg(i)];
%         stnmuse = [stnmuse; stnm(i, :)];
%     end
% end

%% 2. renew vel. model
load('s_wave_para_alxa23.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvp4', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
% load('s_wave_para.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
n=4;
vs = [bestvs1; bestvs2; bestvs4; 4.48];
% vs = [1.1; 3.52; 4.2; 4.48];
vp = [bestvp1; bestvp2; bestvp4; 8.12];
% vp = [bestvp1; bestvp2; 6.75; 8.12];
tk = [besttk1; besttk2+besttk3; besttk4; 0.0];
% dep = bestdep;
dep = 13;
% vs = [1.7; 3.49; 4.3; 4.48];
% vp = [3.3; 5.6; 7.1; 8.12];
% tk= [3.6; 29.2; 16; 0.0];
% dep = 18.1;

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
% pnmod = interp1(xpn, tpn, distsamp, 'spline'); 

% for pPn, model
xppnmin = rayp2 * ( 3*tk(1)/etap1 + (dep-tk(1)+2*tk(2))/etap2 + 2*tk(3)/etap3 );
tppnmin = rayp2 * xppnmin +3*tk(1)*etap1 + (dep-tk(1)+2*tk(2))*etap2 + 2*tk(3)*etap3;
dxppn = linspace(0, 400, nsample2)';
xppn = xppnmin +dxppn;
tppn = tppnmin +rayp2 .*dxppn;
% ppnmod = interp1(xppn, tppn, distsamp, 'spline'); 

% Pn, stack
refdist = 326.0027;
vpn = 8.12;
perip = 1.66;
tp0 = 50.54-perip;                         % time of Pn at refdist 
pnstack = (distsamp - refdist)/vpn + tp0;    % PWS pn curve

% pPn, stack
dtp1 = 5.94;                         % difference time between Pn and pPn
ppnstack = pnstack+dtp1;

% % Pb, stack
% vpb = 6.75;
% peripb = 1.02;
% tp1 = 57.75-peripb;                         % time of Pb at refdist 
% pbstack = (distsamp - refdist)/vpb + tp1;    % PWS pn curve

% Pg picks fitting curve 
a1 = 7.974;
b1 = 46.69;
pgfit = a1.*sqrt((distsamp./b1).^2+1);
raypfit = a1/b1.*distsamp./sqrt(distsamp.^2+b1^2);


%% plot all
% plot P wave

figure
plot(dist, totalshift, 'k*', 'MarkerSize', 10);
set(gca, 'XLIM', [40, 410]);
set(gca, 'XTICK', 40: 50: 410);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);
xlabel('震中距  (/km)', 'Fontsize', 18);
ylabel('偏移量  (/s) ', 'Fontsize', 18);

figure
plot(dist, tpshift, 'k*', 'MarkerSize', 10);
set(gca, 'XLIM', [40, 410]);
set(gca, 'XTICK', 40: 50: 410);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);
xlabel('震中距  (/km)', 'Fontsize', 18);
ylabel('偏移量  (/s) ', 'Fontsize', 18);

figure
plot(distuse, shiftuse, 'k*', 'MarkerSize', 10);

figure
plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;                                                    % pg picks
plot(pgo, dist, 'g.', 'MarkerSize', 8); hold on;                                                    % pg picks

figure
plot(pguse, distuse, 'k.', 'MarkerSize', 8); hold on;                                                    % pg picks
plot(pgo, dist, 'g.', 'MarkerSize', 8); hold on;                                                    % pg picks

