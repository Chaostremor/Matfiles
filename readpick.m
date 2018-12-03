% read time mark of data, to fit the picks to a curve
%

clear;
fid1 = fopen('L:\nodecimate\400nofManPick\timemarknofnod.dat') ;      % strcat�����ַ�������
time = textscan(fid1, '%s %f %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
fclose(fid1) ;
stnm = char(time{1});
dist = time{2};
realpg = time{3};
realsg = time{4};

% Pg picks fitting curve 
distsamp = 0: 0.1: 500;
a1 = 6.104;
b1 = 36.1;
% a1 = 6.443;
% b1 = 39.43;
pgfit = a1.*sqrt((distsamp./b1).^2+1);
raypfit = a1/b1.*distsamp./sqrt(distsamp.^2+b1^2);

a2 = 10.53;
b2 = 37.2;
sgfit = a2.*sqrt((distsamp./b2).^2+1);
raysfit = a2/b2.*distsamp./sqrt(distsamp.^2+b2^2);

% %% 1. read data, lowpass 0.5 hz
% datadir = 'J:\̨ʽ��G��\Alxa\nodecimate\3test400\' ;    % ��������Ŀ¼
% % datadir = 'G:\Alxa\nodecimate\3test400_mis_alxa23\' ;    % ��������Ŀ¼
% fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat�����ַ�������
% weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
% %
% fclose(fid1) ;
% stnm = char(weight{1});     % ̨վ���ǵ�һ��
% [sa, sb] = size(stnm);
% % test to get size of data
% i=1;
% filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡ�ļ�·��
% [ttest, datatest, SAChdrtest] = fget_sac(filename) ;
% disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
% dt = SAChdrtest.times.delta ;
% npts = SAChdrtest.data.trcLen ; 
% %
% [sc, sd] = size(ttest);
% t = zeros(sc,sa);          % t��ʱ�䣬data������
% data1 = zeros(sc,sa);   % data1 ���Z��������
% data2 = zeros(sc, sa);  % data2 ���T��������
% dist = zeros(sa,1);       % dist��ͷ���е����о�
% pg = zeros(sa,1);        % pg��ֱ��P����ʱ
% sg = zeros(sa,1);         % sg��ֱ��S����ʱ
% for i = 1:sa
%     filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡz��������·��
%     [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
%     dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
%     pg(i) = SAChdr.times.t2;
%     sg(i) = SAChdr.times.t3;
%     filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % ��ȡt��������·��
%     [~, data2(:,i), ~] = fget_sac(filename) ; 
% end
% 
% figure
% for i = 1:sa
%     plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(pgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % pg fit
% set(gca, 'XLIM', [0, 100]);
% set(gca, 'XTICK', 0: 10: 100);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);

% figure
% for i = 1:sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % sg fit
% set(gca, 'XLIM', [0, 100]);
% set(gca, 'XTICK', 0: 10: 100);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);

%% 2. read data, bandpass 0.05-0.2 hz
datadir = 'J:\̨ʽ��G��\Alxa\realsynalxa39filter0.05to0.2\' ;    % ��������Ŀ¼
% datadir = 'G:\Alxa\nodecimate\3test400_mis_alxa23\' ;    % ��������Ŀ¼
fid1 = fopen(strcat(datadir,'nweight.dat')) ;      % strcat�����ַ�������
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % ̨վ���ǵ�һ��
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡ�ļ�·��
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
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡz��������·��
    [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
    dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % ��ȡt��������·��
    [~, data2(:,i), ~] = fget_sac(filename) ; 
end

figure
for i = 1:sa
    plot(t(:,i), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
plot(pgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % pg fit
set(gca, 'XLIM', [0, 100]);
set(gca, 'XTICK', 0: 10: 100);
set(gca, 'YLIM', [0, 420]);
set(gca, 'YTick', 0: 20: 420);

% figure
% for i = 1:sa
%     plot(t(:,i), data2(:,i)./(max(abs(data2(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
% end
% plot(sgfit, distsamp, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;            % sg fit
% set(gca, 'XLIM', [0, 150]);
% set(gca, 'XTICK', 0: 10: 150);
% set(gca, 'YLIM', [0, 420]);
% set(gca, 'YTick', 0: 20: 420);