% use data and syn processed by CAP and redo cross-correlation to get sg
% time in real data
%
% Author:
%     C. Song, 2017.8.5
% 
% Version
%     V1.5, use data after alxa35_13
%     <=250 km, cut window
%     >250 km, whole part

clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\test3_13\alxa35_13_' ;    % 数据所在目录
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


% model para.
vs = [2.09; 3.56; 4.12; 4.48];
vp = [3.7; 5.94; 6.30; 8.12];
tk= [6.1; 25.2; 14.9; 0.0];
dep = 13;

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

% cut window for cross-correlation
synstart = zeros(sa, 1);
synend = zeros(sa, 1);
realstart = zeros(sa, 1);
realend = zeros(sa, 1);
synnpt = zeros(sa, 1);
realnpt = zeros(sa, 1);
sgshift = zeros(sa ,1);
flag = 0;


% %%%%%%%   way 2 to cut window
% for i = 1:sa
% %     i = 2;
%     [~, ind2] = min(abs((t2(1: npts(i), i)-synsg(i))));
%     [~, ind1] = min(abs((t1(1: npts(i), i)-synsg(i))));
%     if dist(i) <= 200
%         bdwin = 5;
%         synstart(i) = 1;
%         realstart(i) = 1+ind1-ind2;
%     else
%         ftwin = 8;
%         bdwin = 7;
%         synstart(i) = ind2-round(ftwin./dt);
%         realstart(i) = ind1-round(ftwin./dt);
%     end                
%     synend(i) = ind2+round(bdwin./dt);
% %     npts(i) = ind2+round(win./dt);
%     realend(i) = ind1+round(bdwin./dt);
%     sgshift(i) = t1(ind1, i) - t2(ind2, i);
%     if realstart(i) <=0
%         synstart(i) = 2-realstart(i);
%         realstart(i) = 1;        
%     end
%     synnpt(i) = synend(i) - synstart(i)+1;
%     realnpt(i) = realend(i) - realstart(i)+1;
%     if synnpt(i) ~= realnpt(i)
%         flag = flag+1;
%     end
% end

%%%%%%%%   way 3 to cut window
for i = 1:sa
%     i = 2;
    [~, ind2] = min(abs((t2(1: npts(i), i)-synsg(i))));
    [~, ind1] = min(abs((t1(1: npts(i), i)-synsg(i))));
    if dist(i) <= 250
        bdwin = 5;
        synstart(i) = 1;
        realstart(i) = 1+ind1-ind2;
        synend(i) = ind2+round(bdwin./dt);
        realend(i) = ind1+round(bdwin./dt);
    else       
        synstart(i) = 1;
        realstart(i) = 1;
        synend(i) = npts(i);
        realend(i) = npts(i);
    end                

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

t1use = zeros(800, sa);
t2use = zeros(800, sa);
data1use = zeros(800, sa);
data2use = zeros(800, sa);
for i =1:sa
    t1use(1: nptuse(i), i) = t1(realstart(i): realend(i), i);
    data1use(1: nptuse(i), i) = data1(realstart(i): realend(i), i);
    t2use(1: nptuse(i), i) = t2(synstart(i): synend(i), i);
    data2use(1: nptuse(i), i) = data2(synstart(i): synend(i), i);
    
end

% for i = 1: sa
%     figure (1)
% % i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
% plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
%     plot(synsg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (1)
% end

% i=sa;
% figure
% plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
% plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
% % 
% % figure
% % plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% % plot(t2use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
% % plot(synsg(i), 0, 'b.','MarkerSize', 40);

t3use = zeros(800, sa);
t5use = zeros(800, sa);
maxcoef = zeros(sa, 1);
index = zeros(sa, 1);
nlag = zeros(sa, 1);
tlag = zeros(sa, 1);
realsg = zeros(sa, 1);
for i = 1:sa
%     i=60;
    [coef, lag] = xcorr(data1use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'coeff');
    maxcoef(i) = max(coef);
    index(i) = find(coef==maxcoef(i));
    nlag(i) = lag(index(i));
    tlag(i) = nlag(i)*dt;
    t3use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i);
    t5use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i) + sgshift(i);     % t5 maybe the best !
    realsg(i) =  synsg(i) + tlag(i) + sgshift(i);
end

% i=150;
% figure
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% plot(t2use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
% plot(synsg(i), 0, 'b.','MarkerSize', 40);
% 
% figure
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% plot(t3use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
% plot(t4use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'b'); hold on             % syn data after cut
% plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'm'); hold on             % syn data after cut

% realpg = synpg+tlag;
% realsg = synsg+tlag;

% for i = 1: sa
%     figure (2)
% % i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
%     plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'k'); hold on             % real data after cut
%     plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'r'); hold on             % syn data after cut
%     plot(realsg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (2)
% end

totallag = tlag+sgshift;
figure
plot(dist, totallag, 'b.', 'MarkerSize', 20);
% set(gca, 'YLIM', [-8, 8]);
% set(gca, 'YTick', -8: 2: 8);


%%%%%%%%%% 1. use own crscrl to do pg correction 
distuse = [];
laguse = [];
sguse = [];
% sguse = [];
stnmuse = [];
for i = 1:sa
%     if totallag(i) > -6 && totallag(i) < 0.5
     if dist(i)<=250 && totallag(i) > -5 && totallag(i) < 0.5   
        distuse =  [distuse; dist(i)];
        laguse = [laguse; totallag(i)];
        sguse = [sguse; realsg(i)];
%         sguse = [sguse; realsg(i)];
        stnmuse = [stnmuse; stnm(i, :)];
    end
end

figure
% plot(dist(ind), totalshift(ind), 'b.', 'MarkerSize', 20);
plot(distuse, laguse, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);
