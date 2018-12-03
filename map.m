clear;close all   
%textread函数读取文件，十分强大，好好学习
%读取响应路径的台站信息文件，并去掉若干行“头段”
[lat1,lon1] = textread('G:\qssp2010-code+input\US_Stationlist.txt','%f %f %*[^\n]','headerlines',1);
load coast;
figure (1);hold on;
%axesm lambert;
%等角投影
axesm ('eqdazim','Frame','on','Origin',[38.2963,142.498]);
plotm(lat,long);
scatterm(lat1,lon1,'r','O');
scatterm(38.2963,142.498,'k','^');