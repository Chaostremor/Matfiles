% testinterpolate
%
clear ; clc ; close all;
%% read data, real or synthetic
%cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
%datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录

fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接，filelist事先编辑好
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});
[sa, sb] = size(stnm);
% test to get size of data
i=1;
%filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取文件路径
filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);
data = zeros(sc,sa);
dist = zeros(sa,1);
for i = 1:sa
    %filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取文件路径
    filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end
% plot original data
figure (1)
for i = 1:sa
   %data(:,i) = data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i));
   % normalize
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
% set(gca,'XLim',[0 180]);  % 设置x轴的间隔，长度
% set(gca,'XTick',0:10:180);
% set(gca,'YLim',[30 410]);  % 设置x轴的间隔，长度
% set(gca,'YTick',30:10:410);

for i=1:sa
    if dist(i) > 300
        index = i;
        break
    end
end
nsum = sa-index+1;
tshift = zeros(nsum,1);
%nshift = zeros(nsum,1);
for i = 1:nsum
    tshift(i) = t(1,index)-t(1,i+index-1)+dist(i+index-1)./4.756;     % save two digits of decimal
    shift(i) = tshift(i)/dt;
    
end
dataout=[];
for i=1:nsum
    dout = specshift(data(:,i+index-1),shift(i)-shift(1));     % shift according to the first trace
    dataout = [dataout dout];      % all data has the same length as original
end
tout(:,1) = t(:,index)-tshift(1);      % all data has the same time trace, i.e. first time trace - tshift   
figure (5)
for i = 1:nsum
    plot(dataout(:,i)./(max(abs(dataout(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end