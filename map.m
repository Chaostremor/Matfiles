clear;close all   
%textread������ȡ�ļ���ʮ��ǿ�󣬺ú�ѧϰ
%��ȡ��Ӧ·����̨վ��Ϣ�ļ�����ȥ�������С�ͷ�Ρ�
[lat1,lon1] = textread('G:\qssp2010-code+input\US_Stationlist.txt','%f %f %*[^\n]','headerlines',1);
load coast;
figure (1);hold on;
%axesm lambert;
%�Ƚ�ͶӰ
axesm ('eqdazim','Frame','on','Origin',[38.2963,142.498]);
plotm(lat,long);
scatterm(lat1,lon1,'r','O');
scatterm(38.2963,142.498,'k','^');