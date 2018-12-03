clear ; clc ;
close all;
%% read data, real or synthetic

datadir = 'G:\Alxa\400_15km_real\' ;   % ��������Ŀ¼
fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat�����ַ�������
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊweight{1}  
fclose(fid1) ;
stnm = char(weight{1});
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.r')) ;  % ��ȡ�ļ�·��
%filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % ��ȡ�ļ�·��
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
% load data
[sc, sd] = size(ttest);
treal = zeros(sc,sa);
real = zeros(sc,sa);
dist = zeros(sa,1);
azi = zeros(sa,1);
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.r')) ;  % ��ȡ�ļ�·��
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % ��ȡ�ļ�·��
    [treal(:,i), real(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data�����ݣ�SAChdr��ͷ�α���
    dist(i) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
    azi(i) = SAChdr.evsta.az ;
end

datadir = 'G:\Alxa\400_15km_syn\' ;    % ��������Ŀ¼
%tsyn = zeros(sc,sa);
%syn = zeros(sc,sa);
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.r')) ;  % ��ȡ�ļ�·��
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % ��ȡ�ļ�·��
    [tsyn(:,i), syn(:,i), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data�����ݣ�SAChdr��ͷ�α���
end

fid1 = fopen('G:\Alxa\tshift.txt') ;      % strcat�����ַ�������
tshift = textscan(fid1, '%s %f %f %f %f %f %f %f') ;  % �õ�����һ��cell array,�������Ϊweight{1}  
fclose(fid1) ;
pzshift = tshift{4};
prshift = tshift{5};
for i = 1:sa
    tssyn(:,i) = tsyn(:,i)+prshift(i);
end
figure (1)
for i = 1:sa
    plot(treal(:,i), real(:,i)./(max(abs(real(:,i)))) .*15.0 + double(dist(i)), 'black');hold on;
    plot(tsyn(:,i), syn(:,i)./(max(abs(syn(:,i)))) .*15.0 + double(dist(i)), 'red');hold on;
    plot(tssyn(:,i), syn(:,i)./(max(abs(syn(:,i)))) .*15.0 + double(dist(i)), 'green');hold on;
end
set(gca,'YLIM', [round(min(dist))-10, round(max(dist))+10]);
set(gca,'YLIM', [-10, 500]);