% ��ȡsac�ļ�������������ͷ�α�����ֵɸѡ̨վ�����¿��õ�̨վ����
% ���������ɸѡ�����ǣ����о�distҪ�� [30,90]��֮��
%
clear ; clc ;
cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
datadir = 'G:\newzld\newzld_hinet\' ;   % ��������Ŀ¼
%mkdir('G:\newzld\newzld_hinet_sort\') ;
sortdir = 'G:\newzld\newzld_hinet_sort\' ;   % ɸѡ���Ŀ¼
fid = fopen(strcat(datadir,'filelist')) ;      % strcat�����ַ������ӣ�filelist���ȱ༭��
filelist = textscan(fid, '%s \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
%
fclose(fid) ;
flist = char(filelist{1}) ;   % �ļ����ַ�������  
[sa,sb] = size(flist) ;    % �ļ����ַ�������Ĵ�С
evlat = 42.7568 ;   % �¼�γ��
evlon = 173.0766 ;   % �¼�����
%
j=1 ;
stlon = zeros(sa,1) ;
stlat =  zeros(sa,1) ;
disdeg = zeros(sa,1) ;
azideg = zeros(sa,1) ;
for i = 1:sa
    filename = strcat(datadir,flist(i,:)) ;  % ��ȡ�ļ�·��
    [t, data, SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data�����ݣ�SAChdr��ͷ�α���
    stlon(i) = SAChdr.station.stlo ;   % ע��ͷ�α����Ĵ洢��ʽ
    stlat(i) = SAChdr.station.stla ;
    [disdeg(i), azideg(i)] = distance (evlat, evlon, stlat(i), stlon(i)) ;  % �����㷨��΢С�Ĳ���
    %[diskm, disdeg, bazideg, azideg] = distaz (stlat, stlon, evlat, evlon) ;
    % �ж�����
    if (disdeg(i) >= 30.0) && (disdeg(i) <= 90.0)
        fremlist(j,:) = flist(i,:) ;
        disremin(j,1) = disdeg(i) ;
        copyfile(filename, strcat(sortdir, fremlist(j,:))) ;   % �����ļ���ɸѡ���Ŀ¼
        j=j+1;
    end
end