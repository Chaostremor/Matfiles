% plot travel time fitting inversion result
%
% 2017.6.24, C. Song

clear;
%% plot Pg, Sg arrival time
datadir = 'G:\Alxa\nodecimate\3test400\' ;    % 数据所在目录
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
pg = zeros(sa,1);        % pg是直达P波到时
sg = zeros(sa,1);         % sg是直达S波到时
stla = zeros(sa,1);       % 台站纬度
stlo = zeros(sa,1);       % 台站经度
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取z分量数据路径
    [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data1是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    pg(i) = SAChdr.times.t2;
    sg(i) = SAChdr.times.t3;
    stla(i) = SAChdr.station.stla;
    stlo(i) = SAChdr.station.stlo;
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取t分量数据路径
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end

fid2 = fopen('G:\Alxa\nodecimate\3test400\timemark2.dat') ;      % strcat用于字符串连接
time = textscan(fid2, '%s %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid2) ;
pgo = time{2};     % 台站名是第一列
sgo = time{3};

fid3 = fopen('G:\Alxa\nodecimate\3test400\timemark3.dat') ;      % strcat用于字符串连接
time = textscan(fid3, '%s %f %f %f %f %f\n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid3) ;
shift1 = time{4};     % 台站名是第一列
pg1 = time{5};
sg1 = time{6};

fid4 = fopen('G:\Alxa\nodecimate\3test400\timemark4.dat') ;      % strcat用于字符串连接
time = textscan(fid4, '%s %f %f %f %f %f\n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
fclose(fid4) ;
shift2 = time{4};     % 台站名是第一列
pg2 = time{5};
sg2 = time{6};


distuse = [];
shift1use = [];
pg1use = [];
sg1use = [];
stnmuse = [];
for i = 1:sa
    if shift1(i) > -5 && shift1(i) < 2
        distuse =  [distuse; dist(i)];
        shift1use = [shift1use; shift1(i)];
        pg1use = [pg1use; pg1(i)];
        sg1use = [sg1use; sg1(i)];
        stnmuse = [stnmuse; stnm(i, :)];
    end
end

% figure
% distsamp = 0: 0.1: 500;
% a1 = 7.974;
% b1 = 46.69;
% pgfit = a1.*sqrt((distsamp./b1).^2+1);
% plot(pgo, dist, 'k.', 'MarkerSize', 8); hold on;       % pgo pick
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;            % pg fit
% a1 = 5.32;
% b1 = 32.42;
% pgfit = a1.*sqrt((distsamp./b1).^2+1);
% plot(pg1, dist, 'g.', 'MarkerSize', 8); hold on;       % pgo pick after shift 
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit
% 
% a2 = 10.78;
% b2 = 37.42;
% sgfit = a2.*sqrt((distsamp./b2).^2+1);
% plot(sgo, dist, 'k.', 'MarkerSize', 8); hold on;        % sgo pick
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;            % sg fit
% a2 = 10.16;
% b2 = 36.5;
% sgfit = a2.*sqrt((distsamp./b2).^2+1);
% plot(sg1, dist, 'g.', 'MarkerSize', 8); hold on;        % sgo pick after shift
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % sg fit
% 
% set(gca, 'XLIM', [0, 130]);
% set(gca, 'XTICK', 0: 20: 130);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 50: 420);
% legend('拾取的Pg/Sg到时', '拾取到时拟合曲线', 'cap校正后的到时', 'cap校正后拟合曲线', 'Location', 'Best');
% % legend('拾取的Pg/Sg到时', '拾取到时拟合曲线',  'Location', 'Best');
% xlabel('时间  (/s) ', 'Fontsize', 15);
% ylabel('震中距  (/km)', 'Fontsize', 15);



% %% plot model determined curve on waveform
% figure
% 
% subplot(2, 3, 1);
% for i = 2: 2: sa
%     plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
% end
% load('s_wave_para_alxa23.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvp4', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
% vs = [bestvs1; bestvs2; bestvs4; 4.48];
% vp = [bestvp1; bestvp2; bestvp4; 8.12];
% tk = [besttk1; besttk2+besttk3; besttk4; 0.0];
% dep = bestdep;
% 
% % for direct P, Pg, model 
% nsample1 = 10000;
% rayp1 = linspace(0, 1.0/vp(2)-0.0000005, nsample1)';
% xpg = zeros(nsample1,1);
% tpg = zeros(nsample1,1);
% for ii = 1: nsample1
%     etap1 = sqrt((1.0/vp(1))^2-rayp1(ii)^2);
%     etap2 = sqrt((1.0/vp(2))^2-rayp1(ii)^2);
%     xpg(ii) = rayp1(ii) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
%     tpg(ii) = rayp1(ii)*xpg(ii) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
% end
% 
% % for Pn, model
% nsample2 = 1000;
% rayp2 =  1/vp(4);
% etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
% etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
% etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
% xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
% tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
% dxpn = linspace(0, 400, nsample2)';
% xpn = xpnmin +dxpn;
% tpn = tpnmin +rayp2 .*dxpn;
% 
% % for pPn, model
% xppnmin = rayp2 * ( 3*tk(1)/etap1 + (dep-tk(1)+2*tk(2))/etap2 + 2*tk(3)/etap3 );
% tppnmin = rayp2 * xppnmin +3*tk(1)*etap1 + (dep-tk(1)+2*tk(2))*etap2 + 2*tk(3)*etap3;
% dxppn = linspace(0, 400, nsample2)';
% xppn = xppnmin +dxppn;
% tppn = tppnmin +rayp2 .*dxppn;
% 
% plot(tpg, xpg, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                      % pg model
% text(73, 408, 'Pg', 'fontsize', 10, 'color', 'b');
% plot(tpn, xpn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % pn model
% text(53, 408, 'Pn', 'fontsize', 10, 'color', 'r');
% plot(tppn, xppn, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;                  % ppn model
% text(59, 408, 'pPn', 'fontsize', 10, 'color', 'g');
% set(gca, 'XLIM', [0, 80]);
% set(gca, 'XTICK', 0: 10: 80);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 50: 420);
% xlabel('时间  (/s) ', 'Fontsize', 18);
% ylabel('震中距  (/km)', 'Fontsize', 18);
% 
% subplot(2, 3, 2);
% for i = 2: 2: sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
% end
% % for direct S, Sg
% rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
% xsg = zeros(nsample1,1);
% tsg = zeros(nsample1,1);
% for ii = 1: nsample1
%     etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
%     etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
%     xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
%     tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
% end
% 
% % for Sn
% rays2 = 1/vs(4);
% etas1 = sqrt((1.0/vs(1))^2-rays2^2);
% etas2 = sqrt((1.0/vs(2))^2-rays2^2);
% etas3 = sqrt((1.0/vs(3))^2-rays2^2);
% xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
% tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
% dxsn = linspace(0, 400, nsample2)';
% xsn = xsnmin +dxsn;
% tsn = tsnmin +rays2 .*dxsn;           
% 
% plot(tsg, xsg, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                      % sg model
% text(108, 408, 'Sg', 'fontsize', 10, 'color', 'b');
% plot(tsn, xsn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % sn model
% text(95, 408, 'Sn', 'fontsize', 10, 'color', 'r');
% set(gca, 'XLIM', [0, 130]);
% set(gca, 'XTICK', 0: 10: 130);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 50: 420);
% xlabel('时间  (/s) ', 'Fontsize', 18);
% ylabel('震中距  (/km)', 'Fontsize', 18);
% 
% subplot(2, 3, 3);
% % for i=1: sa
% %     if shift1(i) >= 0
% %         plot(dist(i), shift1(i), 'r*', 'MarkerSize', 10 ); hold on
% %     else
%         plot(dist, shift1, 'k*', 'MarkerSize', 10 ); hold on
% %     end
% % end
% set(gca, 'XLIM', [40, 410]);
% set(gca, 'XTICK', 40: 50: 410);
% set(gca, 'YLIM', [-8, 8]);
% set(gca, 'YTick', -8: 2: 8);
% xlabel('震中距  (/km)', 'Fontsize', 18);
% ylabel('偏移量  (/s) ', 'Fontsize', 18);
% 
% 
% 
% subplot(2, 3, 4);
% for i = 2: 2: sa
%     plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
% end
% load('s_wave_para_alxa28.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvp4', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
% n=4;
% vs = [bestvs1; bestvs2; bestvs4; 4.48];
% vp = [bestvp1; bestvp2; bestvp4; 8.12];
% tk = [besttk1; besttk2+besttk3; besttk4; 0.0];
% dep = bestdep;
% 
% % for direct P, Pg, model 
% nsample1 = 10000;
% rayp1 = linspace(0, 1.0/vp(2)-0.0000005, nsample1)';
% xpg = zeros(nsample1,1);
% tpg = zeros(nsample1,1);
% for ii = 1: nsample1
%     etap1 = sqrt((1.0/vp(1))^2-rayp1(ii)^2);
%     etap2 = sqrt((1.0/vp(2))^2-rayp1(ii)^2);
%     xpg(ii) = rayp1(ii) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
%     tpg(ii) = rayp1(ii)*xpg(ii) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
% end
% 
% % for Pn, model
% nsample2 = 1000;
% rayp2 =  1/vp(4);
% etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
% etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
% etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
% xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
% tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
% dxpn = linspace(0, 400, nsample2)';
% xpn = xpnmin +dxpn;
% tpn = tpnmin +rayp2 .*dxpn;
% 
% % for pPn, model
% xppnmin = rayp2 * ( 3*tk(1)/etap1 + (dep-tk(1)+2*tk(2))/etap2 + 2*tk(3)/etap3 );
% tppnmin = rayp2 * xppnmin +3*tk(1)*etap1 + (dep-tk(1)+2*tk(2))*etap2 + 2*tk(3)*etap3;
% dxppn = linspace(0, 400, nsample2)';
% xppn = xppnmin +dxppn;
% tppn = tppnmin +rayp2 .*dxppn;
% 
% plot(tpg, xpg, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                      % pg model
% text(70, 408, 'Pg', 'fontsize', 10, 'color', 'b');
% plot(tpn, xpn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % pn model
% text(53, 408, 'Pn', 'fontsize', 10, 'color', 'r');
% plot(tppn, xppn, 'linestyle', '-', 'color', 'g', 'LineWidth', 2); hold on;                  % ppn model
% text(59, 408, 'pPn', 'fontsize', 10, 'color', 'g');
% set(gca, 'XLIM', [0, 80]);
% set(gca, 'XTICK', 0: 10: 80);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 50: 420);
% xlabel('时间  (/s) ', 'Fontsize', 18);
% ylabel('震中距  (/km)', 'Fontsize', 18);
% 
% subplot(2, 3, 5);
% for i = 2: 2: sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
% end
% % for direct S, Sg
% rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
% xsg = zeros(nsample1,1);
% tsg = zeros(nsample1,1);
% for ii = 1: nsample1
%     etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
%     etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
%     xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
%     tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
% end
% 
% % for Sn
% rays2 = 1/vs(4);
% etas1 = sqrt((1.0/vs(1))^2-rays2^2);
% etas2 = sqrt((1.0/vs(2))^2-rays2^2);
% etas3 = sqrt((1.0/vs(3))^2-rays2^2);
% xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
% tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
% dxsn = linspace(0, 400, nsample2)';
% xsn = xsnmin +dxsn;
% tsn = tsnmin +rays2 .*dxsn;           
% 
% plot(tsg, xsg, 'linestyle', '-', 'color', 'b', 'LineWidth', 2); hold on;                      % sg model
% text(106, 408, 'Sg', 'fontsize', 10, 'color', 'b');
% plot(tsn, xsn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % sn model
% text(95, 408, 'Sn', 'fontsize', 10, 'color', 'r');
% set(gca, 'XLIM', [0, 130]);
% set(gca, 'XTICK', 0: 10: 130);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 50: 420);
% xlabel('时间  (/s) ', 'Fontsize', 18);
% ylabel('震中距  (/km)', 'Fontsize', 18);
% 
% subplot(2, 3, 6);
% % for i=1: sa
% %     if shift2(i) >= 0
% %         plot(dist(i), shift2(i), 'r*', 'MarkerSize', 10 ); hold on
% %     else
%         plot(dist, shift2, 'k*', 'MarkerSize', 10 ); hold on
% %     end
% % end
% set(gca, 'XLIM', [40, 410]);
% set(gca, 'XTICK', 40: 50: 410);
% set(gca, 'YLIM', [-8, 8]);
% set(gca, 'YTick', -8: 2: 8);
% xlabel('震中距  (/km)', 'Fontsize', 18);
% ylabel('偏移量  (/s) ', 'Fontsize', 18);

pshift = [];
pstnm = [];
plat = [];
plon = [];
nshift = [];
nstnm = [];
nlat = [];
nlon = [];
for i=1: sa
    if shift2(i) >= 0
        pshift = [pshift; shift2(i)];
        pstnm = [pstnm; stnm(i, :)];
        plon = [plon; stlo(i)];
        plat = [plat; stla(i)];
    else
        nshift = [nshift; shift2(i)];
        nstnm = [nstnm; stnm(i, :)];
        nlon = [nlon; stlo(i)];
        nlat = [nlat; stla(i)];
    end
end

























