% earthquake relocation using P travel time
% 
% Author : 
%     C Song, 2017.7.11
% 
%

clear;
% close all;

% datadir = 'G:\Alxa\nodecimate\3test400\' ;    % ��������Ŀ¼
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
% stla = zeros(sa,1);       % ̨վγ��
% stlo = zeros(sa,1);       % ̨վ����
% az = zeros(sa ,1);        % ̨վ��λ��
% gcarc = zeros(sa ,1);   % ̨վ��Բ������
% for i = 1:sa
%     filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % ��ȡz��������·��
%     [t(:,i), data1(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data1�����ݣ�SAChdr��ͷ�α���
%     dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
%     pg(i) = SAChdr.times.t2;
%     sg(i) = SAChdr.times.t3;
%     stla(i) = SAChdr.station.stla;
%     stlo(i) = SAChdr.station.stlo;
%     az(i) = SAChdr.evsta.az;
%     gcarc(i) = SAChdr.evsta.gcarc;
%     filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % ��ȡt��������·��
%     [~, data2(:,i), ~] = fget_sac(filename) ; 
% end
% evla = SAChdr.event.evla;
% evlo = SAChdr.event.evlo;
% save('400km_z_t_data.mat');

load('400km_z_t_data.mat');

maxdist = 200;
for ii=1:sa
    if dist(ii) > maxdist        
        index = ii-1;
        break
    end
end

% check the tlag distribution in 47_18, delete outlier station 15706, 64050, which index are 1 and 12 in nweight.dat 
distuse = cat(1, dist(2: 11), dist(13: index));
stnmuse = cat(1, stnm(2: 11, :), stnm(13: index, :));
stlause = cat(1, stla(2: 11), stla(13: index));
stlouse = cat(1, stlo(2: 11), stlo(13: index));
azuse = cat(1, az(2: 11), az(13: index));
gcarcuse = cat(1, gcarc(2: 11), gcarc(13: index));
figure;
plot(evlo, evla, 'ro'); hold on;
plot(stlouse, stlause, 'k*');

load('s_wave_para_alxa47.mat', 'bestvp1','bestvp2', 'besttk1', 'besttk2');
v1 = bestvp1;
v2 = bestvp2;
h1 = besttk1;
% h2 = besttk2;
h2 = 18-h1;
c = (v1/v2)^2;

size = length(distuse); 
rayp = zeros(size, 1);
G = [];
for ii=1: size
    x = distuse(ii);
    coef(1) = x^2*c;
    coef(2) = -2*x*c*h1;
    coef(3) = x^2*c+c*h1^2-x^2-h2^2;
    coef(4) = -2*x*c*h1+2*x*h1;
    coef(5) = c*h1^2-h1^2;
    b = roots(coef);
    ind = find(imag(b)==0 & real(b)>0);
    trueb = b(ind);
    a = trueb.^2+1;
    rayp(ii) = 1/v1./sqrt(a);
    Gi(1, 1) = rayp(ii)*(evlo-stlouse(ii))*111*cos(deg2rad(evla))/x;
    Gi(1, 2) = rayp(ii)*(evla-stlause(ii))*111/x;
    G = [G; Gi];
end

% fid4 = fopen('G:\Alxa\nodecimate\3test400\timemark4.dat') ;      % strcat�����ַ�������
% time = textscan(fid4, '%s %f %f %f %f %f\n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
% fclose(fid4) ;
% shift2 = time{4};     % ̨վ���ǵ�һ��
% shiftuse = shift2(1: index);
% D = shiftuse;

load('alxa47_18_pg_tlag.mat');
D = cat(1, totallag(2: 11), totallag(13: index));

% m = [dx0; dy0];
M =  (G'*G)\G'*D;
dx0 = M(1,1);
dy0 = M(2,1);
devlo = dx0/(111*cos(deg2rad(evla)));   % unit: degree 
devla = dy0/111;
evlo_reloc = evlo+devlo;
evla_reloc = evla+devla;

figure
plot(distuse, D, 'b.', 'MarkerSize', 20);
set(gca, 'YLIM', [-8, 8]);
set(gca, 'YTick', -8: 2: 8);

save('stat_for_relocation.mat', 'stnmuse');





























