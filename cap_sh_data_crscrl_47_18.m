% use data and syn processed by CAP and redo cross-correlation to get sg
% time in real data
%
% Author:
%     C. Song, 2017.9.5
% 
% Version
%     use data after alxa47_18
%     <=200 km, cut window
%     man s1 (fit)--> alxa45 (cap)--> (45_27 redo crscrl correct)--> s2 (fit)--> alxa47
%     (cap)--> (47_18 redo crscrl correct)--> s3 (fit)


clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\alxa47_18\alxa47_18_' ;    % 数据所在目录
t1 = zeros(800, sa);
t2 = zeros(800, sa);
data1 = zeros(800, sa);
data2 = zeros(800, sa);
npts = zeros(sa ,1);
for i =1:sa
% i=sa;
filename = strcat(datadir,strcat(stnm(i, :), '.0'));   % sh, data
[t, data, ~] = fget_sac(filename) ;
npts(i) = length(data);
t1(1: npts(i), i) = t; 
data1(1: npts(i), i) = data;
filename = strcat(datadir,strcat(stnm(i, :), '.1'));   % sh, syn
[t, data, ~] = fget_sac(filename) ;
t2(1: npts(i), i) = t; 
data2(1: npts(i), i) = data;
end
dt = t2(2, 1) - t2(1, 1); 

% model para.  -->  alxa47_18
vs = [3.28; 3.55; 4.18; 4.48];
vp = [5.50; 5.96; 6.60; 8.12];
tk= [6.7; 31.4; 15.6; 0.0];
dep = 18;

% syn Sg time
nsample1 = 10000;
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

% we need the whole part before pg(syn) and few seconds after it .
% cut window for cross-correlation
maxdist = 200;
for ii = 1:sa
    if dist(ii) > maxdist         % sSn appears after 300 km; Pn > 270km; pPn > 450km ; Sn > ? km
        index = ii-1;
        break
    end
end

synstart = zeros(index, 1);
synend = zeros(index, 1);
realstart = zeros(index, 1);
realend = zeros(index, 1);
synnpt = zeros(index, 1);
realnpt = zeros(index, 1);
sgshift = zeros(index ,1);
flag = 0;

%%%%%%%%   way 3 to cut window
for i = 1:index
%     i = 2;
    [~, ind2] = min(abs((t2(1: npts(i), i)-synsg(i))));
    [~, ind1] = min(abs((t1(1: npts(i), i)-synsg(i))));
    bdwin = 5;
    synstart(i) = 1;
    realstart(i) = 1+ind1-ind2;
    synend(i) = ind2+round(bdwin./dt);
    realend(i) = ind1+round(bdwin./dt);
              
    sgshift(i) = t1(ind1, i) - t2(ind2, i);
    if realstart(i) <=0
        synstart(i) = 2-realstart(i);
        realstart(i) = 1;        
    end
    synnpt(i) = synend(i) - synstart(i)+1;
    realnpt(i) = realend(i) - realstart(i)+1;
    if synnpt(i) ~= realnpt(i)
        flag = flag+1;
    end
end

if flag == 0
    nptuse = realnpt;
end

t1use = zeros(800, index);
t2use = zeros(800, index);
data1use = zeros(800, index);
data2use = zeros(800, index);
for i = 1:index
    t1use(1: nptuse(i), i) = t1(realstart(i): realend(i), i);
    data1use(1: nptuse(i), i) = data1(realstart(i): realend(i), i);
    t2use(1: nptuse(i), i) = t2(synstart(i): synend(i), i);
    data2use(1: nptuse(i), i) = data2(synstart(i): synend(i), i);    
end


t5use = zeros(800, index);
maxcoef = zeros(index, 1);
ind = zeros(index, 1);
nlag = zeros(index, 1);
tlag = zeros(index, 1);
realsg = zeros(index, 1);
for i = 1:index
%     i=sa;
    [coef, lag] = xcorr(data1use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'coeff');
    maxcoef(i) = max(coef);
    ind(i) = find(coef==maxcoef(i));
    nlag(i) = lag(ind(i));
    tlag(i) = nlag(i)*dt;
    t5use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i) + sgshift(i);     % t5 maybe the best !
    realsg(i) =  synsg(i) + tlag(i) + sgshift(i);
end

totallag = tlag+sgshift;
figure
plot(dist(1: index), totallag, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);
set(gca, 'YTick', -8: 2: 8);
set(gca, 'fontsize', 12);
xlabel('Distance  (km) ', 'Fontsize', 15);
ylabel('Shift  (s) ', 'Fontsize', 15);



% i=1;
% figure
% plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
% plot(synsg(i), 0, 'b.','MarkerSize', 40);
% 
% figure
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% plot(t2use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
% plot(synsg(i), 0, 'b.','MarkerSize', 40);

% for i = 1: index
%     figure (10)
% % i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
%     plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'k'); hold on             % real data after cut
%     plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'r'); hold on             % syn data after cut
%     plot(realsg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (10)
% end

fid2 = fopen('G:\Alxa\400nofManPick\timemarknofnod.dat') ;      % strcat用于字符串连接
time1 = textscan(fid2, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid2) ;
pg1 = time1{3};        % man. picked pg/sg time mark
sg1 = time1{4};

% use crscrl result: total lag time + synsg, to correct parts of the
% manually picked time
sg3 = [realsg; sg1(index+1: end)];

% sg1 and fitting curve
figure
distsamp = 0: 0.1: 500;
a2 = 10.53;      % --> alxa45
b2 = 37.2;
sgfit1 = a2.*sqrt((distsamp./b2).^2+1);
plot(sg1, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;            % sg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

% use sg2 to renew the sg curve

% sg2 and fitting curve
figure
distsamp = 0: 0.1: 500;
a2 = 8.355;      % --> alxa??
b2 = 29.44;
sgfit2 = a2.*sqrt((distsamp./b2).^2+1);
% plot(sg2, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % sg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

% sg3 and fitting curve ??
figure
distsamp = 0: 0.1: 500;
a2 = 4.371;      % --> alxa??
b2 = 15.37;
sgfit3 = a2.*sqrt((distsamp./b2).^2+1);
plot(sg3, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sgfit3, distsamp, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;            % sg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

figure
plot(sgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;            % sg fit
plot(sgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % sg fit
plot(sgfit3, distsamp, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;            % sg fit

save('alxa47_18_sg_tlag.mat', 'totallag', 'realsg', 'synsg');


















