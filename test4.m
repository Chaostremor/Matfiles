% check the relationship between real and synthetic data

clear;
%% 1. read data
datadir = 'G:\Alxa\real_syn\' ;    % ��������Ŀ¼
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat�����ַ�������
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % ̨վ���ǵ�һ��
[sa, sb] = size(stnm);
% test to get size of data
i = 1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡ�ļ�·��
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t1 = zeros(sc,sa);          % t��ʱ�䣬data������
data1 = zeros(sc,sa);   % data1 ���Z��������
dist = zeros(sa,1);       % dist��ͷ���е����о�
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡz��������·��
    [t1(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
    dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
end

i = 1;
filename = strcat(datadir,strcat(stnm(i,:),'_syn.z')) ;  % ��ȡ�ļ�·��
[ttest2, datatest2, SAChdrtest2] = fget_sac(filename) ;
[sc2, ~] = size(ttest2);
t2 = zeros(sc2, sa);
data2 = zeros(sc2, sa);
firstp = zeros(sa,1);        % pg��ֱ��P����ʱ
firsts = zeros(sa,1);         % sg��ֱ��S����ʱ
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.z')) ;  % ��ȡz��������·��
    [t2(:,i), data2(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
    firstp(i) = SAChdr.times.t1;
    firsts(i) = SAChdr.times.t2;
end

fid2 = fopen('G:\Alxa\real_syn\timeshiftfromcap28_18') ;      % strcat�����ַ�������
time = textscan(fid2, '%s %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
%
fclose(fid2) ;
tpshift = time{2};     % ̨վ���ǵ�һ��
corshift = time{3};
totalshift = tpshift+corshift;

m1 = 35;
m2 = 80;
f2_pnl = 0.2;
cuts = firstp-0.3*m1+tpshift;
cute = firsts+tpshift;
for i=1:sa
    if cute(i) < (cuts(i)+2/f2_pnl+0.3*m1)
        cute(i) = cuts(i)+2/f2_pnl+0.3*m1;
    end
end

load('s_wave_para_alxa28.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvp4', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
% load('s_wave_para.mat', 'besttk1', 'besttk2', 'besttk3', 'besttk4', 'bestvp1', 'bestvp2', 'bestvs1', 'bestvs2', 'bestvs4', 'bestdep');
n=4;
vs = [bestvs1; bestvs2; bestvs4; 4.48];
% vs = [1.1; 3.52; 4.2; 4.48];
vp = [bestvp1; bestvp2; bestvp4; 8.12];
% vp = [bestvp1; bestvp2; 6.75; 8.12];
tk = [besttk1; besttk2+besttk3; besttk4; 0.0];
dep = 18;
% dep = 13.1;
% vs = [1.7; 3.49; 4.3; 4.48];
% vp = [3.3; 5.6; 7.1; 8.12];
% tk= [3.6; 29.2; 16; 0.0];
% dep = 18.1;

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
synpgint = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
realpgint = synpgint+totalshift;

figure
% for i = 1:sa
i = 2;
% ./(max(abs(data1(:,i)))) .*15.0
    plot(t1(:, i), data1(:, i) + double(dist(i)), 'LineWidth',2, 'color', 'k');hold on;
    plot(t2(:, i), data2(:, i) + double(dist(i)), 'LineWidth',2, 'color', 'r');hold on;
    plot(synpgint(i), double(dist(i)), 'g.','MarkerSize', 20); hold on;
%     plot(firstp(i), double(dist(i)), 'k.','MarkerSize', 20); hold on;
    plot(cuts(i), double(dist(i)), 'b.', 'MarkerSize', 20); hold on;
    plot(cute(i), double(dist(i)), 'b.', 'MarkerSize', 20); hold on;
% end
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);

t3 = zeros(sc2, sa);
for i=1:sa
    t3(:, i) = t2(:, i)+totalshift(i);
end

figure
% for i = 1:sa
i = 2;
% ./(max(abs(data1(:,i)))) .*15.0
    plot(t1(:, i), data1(:, i) + double(dist(i)), 'LineWidth',2, 'color', 'k');hold on;
    plot(t3(:, i), data2(:, i) + double(dist(i)), 'LineWidth',2, 'color', 'r');hold on;

%     plot(cuts(i), double(dist(i)), 'b.', 'MarkerSize', 20); hold on;
%     plot(cute(i), double(dist(i)), 'b.', 'MarkerSize', 20); hold on;
% end
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);

