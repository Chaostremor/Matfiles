% use data and syn processed by CAP and redo cross-correlation to get pg
% time in real data
%
% Author:
%     C. Song, 2017.9.2
% 
% Version
%     use data after alxa45_27
%     <=200 km, cut window
%     man p1 (fit)--> alxa45 (cap)--> (redo crscrl correct)--> p2 (fit)--> alxa47


clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\alxa45_27\alxa45_27_' ;    % 数据所在目录
t1 = zeros(350, sa);
t2 = zeros(350, sa);
data1 = zeros(350, sa);
data2 = zeros(350, sa);
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

% model para.  -->  alxa45_27
vs = [2.05; 3.64; 4.16; 4.48];
vp = [3.6; 6.08; 6.50; 8.12];
tk= [5.8; 26.8; 16.1; 0.0];
dep = 27;

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
pgshift = zeros(index ,1);
flag = 0;

%%%%%%%%   way 3 to cut window

for i = 1:index
%     i = 2;
    [~, ind2] = min(abs((t2(1: npts(i), i)-synpg(i))));
    [~, ind1] = min(abs((t1(1: npts(i), i)-synpg(i))));
    bdwin = 3.5;
    synstart(i) = 1;
    realstart(i) = 1+ind1-ind2;
    synend(i) = ind2+round(bdwin./dt);
    realend(i) = ind1+round(bdwin./dt);

    pgshift(i) = t1(ind1, i) - t2(ind2, i);
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

t1use = zeros(350, index);
t2use = zeros(350, index);
data1use = zeros(350, index);
data2use = zeros(350, index);
for i = 1:index
    t1use(1: nptuse(i), i) = t1(realstart(i): realend(i), i);
    data1use(1: nptuse(i), i) = data1(realstart(i): realend(i), i);
    t2use(1: nptuse(i), i) = t2(synstart(i): synend(i), i);
    data2use(1: nptuse(i), i) = data2(synstart(i): synend(i), i);    
end


t5use = zeros(350, index);
maxcoef = zeros(index, 1);
ind = zeros(index, 1);
nlag = zeros(index, 1);
tlag = zeros(index, 1);
realpg = zeros(index, 1);
for i = 1:index
%     i=sa;
    [coef, lag] = xcorr(data1use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'coeff');
    maxcoef(i) = max(coef);
    ind(i) = find(coef==maxcoef(i));
    nlag(i) = lag(ind(i));
    tlag(i) = nlag(i)*dt;
    t5use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i) + pgshift(i);     % t5 maybe the best !
    realpg(i) =  synpg(i) + tlag(i) + pgshift(i);
end

totallag = tlag+pgshift;
figure
% ind = find(dist>=250);
% plot(dist(ind), totallag(ind), 'b.', 'MarkerSize', 20);
plot(dist(1: index), totallag, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);
set(gca, 'fontsize', 12);
xlabel('Distance  (km) ', 'Fontsize', 15);
ylabel('Shift  (s) ', 'Fontsize', 15);

% for i = 1: index
%     figure (10)
% % i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
%     plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'k'); hold on             % real data after cut
%     plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'r'); hold on             % syn data after cut
%     plot(realpg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (10)
% end

fid2 = fopen('G:\Alxa\400nofManPick\timemarknofnod.dat') ;      % strcat用于字符串连接
time1 = textscan(fid2, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid2) ;
pg1 = time1{3};        % man. picked pg/sg time mark
sg1 = time1{4};

% use crscrl result: total lag time + synpg, to correct parts of the
% manually picked time
pg2 = [realpg; pg1(index+1: end)];

save('alxa45_27_pg_tlag.mat', 'totallag', 'pg2', 'dist');

figure
distsamp = 0: 0.1: 500;
a1 = 6.104;      % alxa45
b1 = 36.1;
pgfit1 = a1.*sqrt((distsamp./b1).^2+1);
plot(pg1, dist, 'k.', 'MarkerSize', 8); hold on;
plot(pgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

% use pg2 to renew the pg curve

figure
distsamp = 0: 0.1: 500;
a1 = 4.64;      % --> alxa??
b1 = 27.37;
pgfit2 = a1.*sqrt((distsamp./b1).^2+1);
plot(pg2, dist, 'k.', 'MarkerSize', 8); hold on;
plot(pgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

figure
plot(pgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;            % pg fit
plot(pgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit











