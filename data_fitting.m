% fit the Sg, Pg arrival time in least-square sense 
% eg. Sg, Pg
% 
% Author: C. Song,  2017.4.17
% 
clear ; close all;
%% read data, real or synthetic
datadir = 'G:\Alxa\nodecimate\3test400\' ;    % 数据所在目录
%fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接
fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
%filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data = zeros(sc,sa);
dist = zeros(sa,1);       % dist是头段中的震中距
pg = zeros(sa,1);        % pg是直达P波到时
sg = zeros(sa,1);         % sg是直达S波到时
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    pg(i) = SAChdr.times.t2;
    sg(i) = SAChdr.times.t3;
end
figure (1)
for i = 1:sa
   % normalize
   %plot(t(:,i), data(:,i)*100 + double(dist(i)), 'LineWidth',1);hold on;
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
set(gca, 'XLIM', [0, 300]);
set(gca, 'XTICK', 0:20:300);
set(gca, 'YLIM', [round(min(dist))-10, round(max(dist))+10]);
set(gca, 'YTick', round(min(dist))-10:10:round(max(dist))+10);

save('pgsgpicktime.mat', 'pg', 'sg', 'dist');

% ft = fittype('a.*sqrt(1+(x./b).^2)', 'coeff', {'a', 'b'}, 'indep', 'x');
% [pgcurve, goodness] = fit(dist, pg, ft);
% figure (3)
% plot(pgcurve, pg, dist);

% use curve fitting tool
% x data: pg/ sg
% y data: dist
% custom equation: y=f(x)=b.*sqrt((x./a).^2-1)=> dist = b.*sqrt((T./a).^2-1)
figure (2)
%scatter(pg, dist, 8, 'black', 'filled'); hold on;
%scatter(sg, dist, 8, 'black', 'filled'); hold on;
plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sg, dist, 'k.', 'MarkerSize', 8); hold on;
a1 = 5.969;
b1 = 35.05;
vpg = b1/a1;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);
plot(pgcurve, distcurve, 'b-'); hold on;
a2 = 10.34;
b2 = 35.95;
vsg = b2/a2;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
plot(sgcurve, distcurve, 'b-');

figure (3)
for i = 1:sa
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sg, dist, 'k.', 'MarkerSize', 8); hold on;
pgcurve = a1.*sqrt((dist./b1).^2+1);
plot(pgcurve, dist, 'b-', 'LineWidth', 2); hold on;
sgcurve = a2.*sqrt((dist./b2).^2+1);
plot(sgcurve, dist, 'b-', 'LineWidth', 2);
set(gca, 'XLIM', [0, 300]);
set(gca, 'XTICK', 0:20:300);
set(gca, 'YLIM', [round(min(dist))-10, round(max(dist))+10]);
set(gca, 'YTick', round(min(dist))-10:10:round(max(dist))+10);

thick = 0.1:0.1:10;
ntk = length(thick);
h = 10:0.1:30;
nh = length(h);
vp = 3.0:0.01:5.5;
nvp = length(vp);
nsample =10000;
pP = linspace(0, 1.0/5.87-0.00005, nsample)';
hflag = 2;
misfit = zeros(ntk, nh, nvp);
for i = 1:ntk                     % thickness of sediment, thick
    for j = 1:nh                  % source depth, h
        for k = 1:nvp           % vp of sediment, vp
            XP=zeros(nsample,1);
            TP=zeros(nsample,1);
            for m = 1:nsample
                 etaP1 = sqrt((1.0/vp(k))^2-pP(m)^2);
                 etaP2 = sqrt((1.0/5.87)^2-pP(m)^2);
                 XP(m) = pP(m) * (thick(i)/etaP1+(h(j)-thick(i))/etaP2);
                 TP(m) = pP(m)*XP(m)+thick(i)*etaP1+(h(j)-thick(i))*etaP2;
            end
            pgfit = interp1(XP, TP, dist, 'spline');
            misfit(i, j, k) = sum((pgfit-pg).^2);
        end
    end
end
[bestfit, index] = min(misfit(:));
[itk, jh, kvp] = ind2sub([ntk, nh, nvp], index);
besttk = thick(itk);
besth = h(jh);
bestvp = vp(kvp);
XP=zeros(nsample,1);
TP=zeros(nsample,1);
for m = 1:nsample
    etaP1 = sqrt((1.0/bestvp)^2-pP(m)^2);
    etaP2 = sqrt((1.0/5.87)^2-pP(m)^2);
    XP(m) = pP(m) * (besttk/etaP1+(besth-besttk)/etaP2);
    TP(m) = pP(m)*XP(m)+besttk*etaP1+(besth-besttk)*etaP2;
end
pgint = interp1(XP, TP, distcurve, 'spline');

figure(4)
plot(pgcurve, dist, 'b-'); hold on;
plot(TP, XP, 'k-'); hold on;
plot(pgint, distcurve, 'r-'); 




