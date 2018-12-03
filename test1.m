% fit model check
%
% USAGE:
%     check if the result model from travel_time_curve_fitting is reasonable,    
%     plot the travel time curve on the data based on the new model
%
%  Author:   C. Song, 2017.5.18
%       

%% 1. read data
datadir = 'G:\Alxa\synnofilter\' ;    % ��������Ŀ¼
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat�����ַ�������
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % ̨վ���ǵ�һ��
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'_syn.z')) ;  % ��ȡ�ļ�·��
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t��ʱ�䣬data������
data1 = zeros(sc,sa);   % data1 ���Z��������
data2 = zeros(sc, sa);  % data2 ���T��������
dist = zeros(sa,1);       % dist��ͷ���е����о�
p = zeros(sa,1);        % pg��ֱ��P����ʱ
s = zeros(sa,1);         % sg��ֱ��S����ʱ
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.z')) ;  % ��ȡz��������·��
    [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
    dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
    p(i) = SAChdr.times.t1;
    s(i) = SAChdr.times.t2;
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % ��ȡt��������·��
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end

n=4;
vs = [1.1; 3.52; 3.9; 4.48];
vp = [2.5; 5.53; 6.1; 8.12];
tk= [2.3; 22.2; 18.8; 0.0];
dep = 15;

%% 3. calculate P wave phase travel time curve
distsamp = 0: 0.1: 500;

% for direct P, Pg, model 
nsample1 = 10000;
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for ii = 1: nsample1
    etap1 = sqrt((1.0/vp(1))^2-rayp1(ii)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(ii)^2);
    xpg(ii) = rayp1(ii) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
    tpg(ii) = rayp1(ii)*xpg(ii) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
% pgmod = interp1(xpg, tpg, distsamp, 'spline');       % model determined pg curve

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


% Pg picks fitting curve 
a1 = 6.188;
b1 = 33.83;
pgfit = a1.*sqrt((distsamp./b1).^2+1);
raypfit = a1/b1.*distsamp./sqrt(distsamp.^2+b1^2);

%% 4. calculate S wave phase travel time curve

% for direct S, Sg
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for ii = 1: nsample1
    etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
    xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
    tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
% sgmod = interp1(xsg, tsg, distsamp, 'spline');

% for Sn
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
dxsn = linspace(0, 400, nsample2)';
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;           
% snmod = interp1(xsn, tsn, distsamp, 'spline');       % model determined sn curve

% Sn, stack
refdist = 326.0027;
vsn = 4.48;
peris = 1.46;
ts0 = 85.96-peris;                         % time of Sn at refdist 
snstack = (distsamp - refdist)/vsn + ts0;    % PWS sn curve

% Sg picks fitting curve
a2 = 10.01;
b2 = 34.37;
sgfit = a2.*sqrt((distsamp./b2).^2+1);
raysfit = a2/b2.*distsamp./sqrt(distsamp.^2+b2^2);


%% plot all
% plot P wave
figure
for i = 1:sa
%     temp1 = data(:, i)./(max(abs(data(:, i)))) .*30.0 ;
%     temp2 = (0.5*(temp1+abs(temp1)));
%     plot(t(:, i), temp1+ double(dist(i)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
%     fill(t(:, i), temp2+ double(dist(i)), 'r', 'edgealpha', 0); hold on;
    plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
end
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;            % pg fit
% plot(pnstack, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;       % pn stack
% plot(ppnstack, distsamp, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;       % ppn stack
% % plot(pbstack, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;       % pb stack
% plot(tpg, xpg, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % pg model
% plot(tpn, xpn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % pn model
% plot(p, dist, 'g.', 'MarkerSize', 8); hold on;                                                    % pg picks
% plot(tppn, xppn, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                  % ppn model
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 20: 420);

% plot S wave
figure
for i = 1:sa
%     temp1 = data(:, i)./(max(abs(data(:, i)))) .*30.0 ;
%     temp2 = (0.5*(temp1+abs(temp1)));
%     plot(t(:, i), temp1+ double(dist(i)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
%     fill(t(:, i), temp2+ double(dist(i)), 'r', 'edgealpha', 0); hold on;
    plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'linestyle', '-', 'color', 'k', 'LineWidth',1);hold on;
end

% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % sg fit
% plot(snstack, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;       % sn stack
plot(tsg, xsg, 'linestyle', '-', 'color', 'r', 'LineWidth', 2); hold on;                      % sg model
plot(tsn, xsn, 'linestyle', '--', 'color', 'r', 'LineWidth', 2); hold on;                      % sn model
plot(s, dist, 'g.', 'MarkerSize', 8); hold on;                                                    % sg picks
set(gca, 'XLIM', [0, 150]);
set(gca, 'XTICK', 0: 10: 150);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 20: 420);

