% use data and syn processed by CAP and redo cross-correlation to get pg
% time in real data
%
% Author:
%     C. Song, 2017.8.4
% 
% Version
%     V4, use data after alxa37_16 or alxa39_16,
%     manual pick-->alxa35(dt=0.1)-->CAP-->alxa35_13-->crscrl corrected
%     time-->alxa37/39-->CAP-->alxa37_16/39_16

clear;
datadir = 'G:\Alxa\real_syn\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
dist = weight{2};

datadir = 'G:\Alxa\alxa39_16\alxa39_16_' ;    % 数据所在目录
% datadir = 'G:\Alxa\alxa37_16\alxa37_16_' ;    % 数据所在目录
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

fid2 = fopen('G:\Alxa\alxa39_16\timeshiftfromcap39_16') ;      % strcat用于字符串连接
% fid2 = fopen('G:\Alxa\alxa37_16\timeshiftfromcap37_16') ;      % strcat用于字符串连接
capshift = textscan(fid2, '%s %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid2) ;
conshft = capshift{2};     % 台站名是第一列
corshftp = capshift{3};
totalshft = conshft+corshftp;

% model para.
vs = [2.52; 3.65; 4.04; 4.48];
vp = [3.8; 6.18; 6.40; 8.12];
tk= [3.7; 31.6; 17.4; 0.0];
dep = 16;

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

% % syn Sg time
% rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
% xsg = zeros(nsample1,1);
% tsg = zeros(nsample1,1);
% for ii = 1: nsample1
%     etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
%     etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
%     xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
%     tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
% end
% synsg = interp1(xsg, tsg, dist, 'spline');

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

% %%%%%%%% way 1 to cut window
% for i = 1:sa
% %     i = 2;
% %     win = 3.0;                          % proper value
%     if dist(i) >= 0 && dist(i) <180
%         win = 3.0;
%     elseif dist(i) >= 180 && dist(i) <250
%         win = 1.6;
%     else
%         win = 0.4;
%     end
%     [~, ind2] = min(abs((t2(1: npts(i), i)-synpg(i))));
%     synstart(i) = 1;
%     synend(i) = ind2+round(win./dt);
% %     npts(i) = ind2+round(win./dt);
%     [~, ind1] = min(abs((t1(1: npts(i), i)-synpg(i))));
%     realstart(i) = 1+ind1-ind2;
%     realend(i) = ind1+round(win./dt);
%     pgshift(i) = t1(ind1, i) - t2(ind2, i);
%     if realstart(i) <=0
%         synstart(i) = 2-realstart(i);
%         realstart(i) = 1;        
%     end
%     synnpt(i) = synend(i) - synstart(i)+1;
%     realnpt(i) = realend(i) - realstart(i)+1;
%     if synnpt ~= realnpt(i)
%         flag = flag+1;
%     end
% end

%%%%%%%%   way 2 to cut window
for i = 1:sa
%     i = 2;
    [~, ind2] = min(abs((t2(1: npts(i), i)-synpg(i))));
    [~, ind1] = min(abs((t1(1: npts(i), i)-synpg(i))));
    if dist(i) <= 200
        bdwin = 3.5;
        synstart(i) = 1;
        realstart(i) = 1+ind1-ind2;
    else
%         ftwin = 6;
%         bdwin = 2.5;
        ftwin = 7;
        bdwin = 1.5;
        synstart(i) = ind2-round(ftwin./dt);
        realstart(i) = ind1-round(ftwin./dt);
    end                
    synend(i) = ind2+round(bdwin./dt);
%     npts(i) = ind2+round(win./dt);
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

% %%%%%%%%   way 3 to cut window
% for i = 1:sa
% %     i = 2;
%     [~, ind2] = min(abs((t2(1: npts(i), i)-synpg(i))));
%     [~, ind1] = min(abs((t1(1: npts(i), i)-synpg(i))));
%     if dist(i) <= 250
%         bdwin = 3.5;
%         synstart(i) = 1;
%         realstart(i) = 1+ind1-ind2;
%         synend(i) = ind2+round(bdwin./dt);
%         realend(i) = ind1+round(bdwin./dt);
%     else       
%         synstart(i) = 1;
%         realstart(i) = 1;
%         synend(i) = npts(i);
%         realend(i) = npts(i);
%     end                
% 
%     pgshift(i) = t1(ind1, i) - t2(ind2, i);
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


if flag == 0
    nptuse = realnpt;
end

t1use = zeros(175, sa);
t2use = zeros(175, sa);
data1use = zeros(175, sa);
data2use = zeros(175, sa);
for i =1:sa
    t1use(1: nptuse(i), i) = t1(realstart(i): realend(i), i);
    data1use(1: nptuse(i), i) = data1(realstart(i): realend(i), i);
    t2use(1: nptuse(i), i) = t2(synstart(i): synend(i), i);
    data2use(1: nptuse(i), i) = data2(synstart(i): synend(i), i);
    
end

% for i = 76: sa
    figure
i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
%     plot(synpg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (1)
% end

% for i = 76: sa
%     figure (1)
% % i = 2;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% plot(t2use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
%     plot(synpg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (1)
% end

i=sa;
figure
plot(t1(1: npts(i), i), data1(1: npts(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data
plot(t2(1: npts(i), i), data2(1: npts(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data
plot(synpg(i), 0, 'b.', 'MarkerSize', 40);
plot(synpn(i), 0, 'y.', 'MarkerSize', 40);
plot(synppn(i), 0, 'm.', 'MarkerSize', 40);
% 
% figure
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% plot(t2use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut

t3use = zeros(175, sa);
% t4use = zeros(175, sa);
t5use = zeros(175, sa);
maxcoef = zeros(sa, 1);
index = zeros(sa, 1);
nlag = zeros(sa, 1);
tlag = zeros(sa, 1);
realpg = zeros(sa, 1);
for i = 1:sa
%     i=sa;
    [coef, lag] = xcorr(data1use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'coeff');
    maxcoef(i) = max(coef);
    index(i) = find(coef==maxcoef(i));
    nlag(i) = lag(index(i));
    tlag(i) = nlag(i)*dt;
    t3use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i);
%     t4use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i) + tpshift(i);
    t5use(1: nptuse(i), i) = t2use(1: nptuse(i), i) + tlag(i) + pgshift(i);     % t5 maybe the best !
    realpg(i) =  synpg(i) + tlag(i) + pgshift(i);
end

% i=7;
% figure
% plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth',2, 'color', 'k'); hold on             % real data after cut
% % plot(t3use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'r'); hold on             % syn data after cut
% % plot(t4use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'b'); hold on             % syn data after cut
% plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth',2, 'color', 'm'); hold on             % syn data after cut
% plot(realpg(i), 0, 'b.', 'MarkerSize', 40);


% realpg = synpg+tlag;
% realsg = synsg+tlag;

% for i = 50: sa
%     figure (2)
figure
i = sa;
% i
% % ./max(data1(1: npts(i), i))*10, ./max(data2(1: npts(i), i))*10
    plot(t1use(1: nptuse(i), i), data1use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'k'); hold on             % real data after cut
    plot(t5use(1: nptuse(i), i), data2use(1: nptuse(i), i), 'LineWidth', 2, 'color', 'r'); hold on             % syn data after cut
    plot(realpg(i), 0, 'b.', 'MarkerSize', 40);
%     pause (5);
%     close (2)
% end

totallag = tlag+pgshift;
figure
plot(dist, totallag, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);

% figure
% plot(dist, totalshft, 'b.', 'MarkerSize', 20);
% set(gca, 'YLIM', [-8, 8]);
% set(gca, 'YTick', -8: 2: 8);

distuse = [];
laguse = [];
pguse = [];
% sguse = [];
stnmuse = [];
for i = 1:sa
    if totallag(i) > -2.5 && totallag(i) < 1
        distuse =  [distuse; dist(i)];
        laguse = [laguse; totallag(i)];
        pguse = [pguse; realpg(i)];
%         sguse = [sguse; realsg(i)];
        stnmuse = [stnmuse; stnm(i, :)];
    end
end

