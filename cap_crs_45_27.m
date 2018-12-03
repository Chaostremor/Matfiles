% read cap crscrl result from model alxa45, depth 27km

clear;
fid1 = fopen('G:\Alxa\400nofManPick\timemarknofnod.dat') ;      % strcat用于字符串连接
time1 = textscan(fid1, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid1) ;
stnm = char(time1{1});
dist = time1{2};
sa = length(dist);
pg1 = time1{3};        % man. picked pg/sg time mark
sg1 = time1{4};
% 
% fid2 = fopen('G:\Alxa\timeshiftfromcap45_27') ;      % strcat用于字符串连接
% time2 = textscan(fid2, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
% fclose(fid2) ;
% conshft = time2{2};
% corshftp = time2{3};
% corshfts = time2{4};
% shftp = conshft+corshftp;
% shfts = conshft+corshfts;
% 
% figure
% plot(dist, shftp, 'k.', 'MarkerSize', 8); hold on;        % pnl shift between real and syn of alxa45_27
% 
% figure
% plot(dist, shfts, 'r.', 'MarkerSize', 8); hold on;        % sw shift between real and syn of alxa45_27

figure
distsamp = 0: 0.1: 500;
% a1 = 7.974;       % alxa23
% b1 = 46.69; 
% pgfit = a1.*sqrt((distsamp./b1).^2+1);   
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;            % pg fit
% a1 = 5.32;        % alxa28 
% b1 = 32.42;
% pgfit = a1.*sqrt((distsamp./b1).^2+1);
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit
a1 = 6.104;      % alxa45
b1 = 36.1;
pgfit1 = a1.*sqrt((distsamp./b1).^2+1);
plot(pg1, dist, 'k.', 'MarkerSize', 8); hold on;
plot(pgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

fid3 = fopen('G:\Alxa\timeshiftfromcap45_27') ;      % strcat用于字符串连接
time3 = textscan(fid3, '%s %f %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid3) ;
conshft = time3{2};
corshftp = time3{3};
corshfts = time3{4};
shftp = conshft+corshftp;
shfts = conshft+corshfts;

% figure
% plot(dist, conshft, 'k.', 'MarkerSize', 8); hold on;       % constant shift
% 
% figure
% plot(dist, corshftp, 'k.', 'MarkerSize', 8); hold on;      % cross-correlation shift

%% correct pg
figure
plot(dist, shftp, 'k.', 'MarkerSize', 8); hold on;        % pnl shift between real and syn of alxa45_27

% model para.
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

maxdist = 200;
for ii = 1:sa
    if dist(ii) > maxdist         % sSn appears after 300 km; Pn > 270km; pPn > 450km ; Sn > ? km
        index = ii-1;
        break
    end
end
pg2 = zeros(index,1);
for ii = 1:index
    pg2(ii) = synpg(ii)+shftp(ii);
end
pg2 = [pg2; pg1(index+1: end)];

figure
distsamp = 0: 0.1: 500;
a1 = 5.554;      % alxa45
b1 = 32.62;
pgfit2 = a1.*sqrt((distsamp./b1).^2+1);
plot(pg2, dist, 'k.', 'MarkerSize', 8); hold on;
plot(pgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

figure
plot(pgfit1, distsamp, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;
plot(pgfit2, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;
set(gca, 'XLIM', [0, 130]);
set(gca, 'XTICK', 0: 20: 130);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 50: 420);

%% correct sg
figure
plot(dist, shfts, 'k.', 'MarkerSize', 8); hold on;        % sw shift between real and syn of alxa45_27

% syn Sg time
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

% 从上图可以看出，sw部分过长，直接用cap的偏移时间作sg的校正误差过大，不采纳








