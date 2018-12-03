% Relocation used Pg lag time data from final model 47_18, after that we
% get alxa49_18, we need to check whether the lag time is improved or not
% 
%
% Author:
%     C. Song, 2017.9.6
% 
% Version
%     use data after alxa49_18
%     <=200 km, cut window
%     man p1 (fit)--> alxa45 (cap)--> (45_27 redo crscrl correct)--> p2 (fit)--> alxa47
%     (cap)--> (47_18 redo crscrl correct)--> relocation --> reprocess data
%     --> (cap)--> 49_18 redo crscrl correct


clear;
datadir = 'G:\Alxa\alxa49_18\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'reweight_reloc.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\alxa49_18\alxa49_18_' ;    % 数据所在目录
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

% model para.  -->  alxa49_18, model 49 is same as 47
vs = [3.28; 3.55; 4.18; 4.48];
vp = [5.50; 5.96; 6.60; 8.12];
tk= [6.7; 31.4; 15.6; 0.0];
dep = 18;

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

synstart = zeros(sa, 1);
synend = zeros(sa, 1);
realstart = zeros(sa, 1);
realend = zeros(sa, 1);
synnpt = zeros(sa, 1);
realnpt = zeros(sa, 1);
pgshift = zeros(sa ,1);
flag = 0;

%%%%%%%%   way 3 to cut window

for i = 1:sa
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

t1use = zeros(350, sa);
t2use = zeros(350, sa);
data1use = zeros(350, sa);
data2use = zeros(350, sa);
for i = 1:sa
    t1use(1: nptuse(i), i) = t1(realstart(i): realend(i), i);
    data1use(1: nptuse(i), i) = data1(realstart(i): realend(i), i);
    t2use(1: nptuse(i), i) = t2(synstart(i): synend(i), i);
    data2use(1: nptuse(i), i) = data2(synstart(i): synend(i), i);    
end


t5use = zeros(350, sa);
maxcoef = zeros(sa, 1);
ind = zeros(sa, 1);
nlag = zeros(sa, 1);
tlag = zeros(sa, 1);
realpg = zeros(sa, 1);
for i = 1:sa
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
plot(dist(1: sa), totallag, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);

save('alxa49_18_pg_tlag.mat', 'totallag');

% figure
% i = 13;
% % i
% plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
% plot(synpg(i), 0, 'b.', 'MarkerSize', 40);
% 
% figure (10)
% i = 13;
% % i
%     plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'k'); hold on             % real data after cut
%     plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'r'); hold on             % syn data after cut
%     plot(realpg(i), 0, 'b.', 'MarkerSize', 40);





