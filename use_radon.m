% use Radon transform
% author:  C. Song,  2017.3.15
%
clear ; clc ; close all;
%% read data, real or synthetic
%cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
%datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录

fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接，filelist事先编辑好
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取文件路径
%filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
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
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end
% plot original data
figure (1)
for i = 1:sa
   % normalize
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
% set(gca,'XLim',[0 180]);  % 设置x轴的间隔，长度
% set(gca,'XTick',0:10:180);
% set(gca,'YLim',[30 410]);  % 设置x轴的间隔，长度
% set(gca,'YTick',30:10:410);

%% there are 3 ways to shift the data according to one specific p parameter ( or velocity )
%% 1. resample
% newrate = 100;
% for i = 1:sa
%     [datainter(:,i), tinter(:,i)] = resample(data(:,i), t(:,i), newrate, 'spline');
% end
% shift and plot
% for i=1:sa
%     if dist(i) > 300
%         index = i;
%         break
%     end
% end
% [newnpts, ttt] = size(tinter(:,1));       % ttt is useless
% nsum = sa-index+1;
% tshift = zeros(nsum,1);
% nshift = zeros(nsum,1);
% for i = 1:nsum
%     tshift(i) = roundn(dist(i+index-1)./4.756, -2);     % save two digits of decimal
%     nshift(i) = round(tshift(i)/newdt);
% end
% figure (2)
% for i = 1:nsum
%     tinter(:,i+index-1) =  tinter(:,i+index-1) - tshift(i);
%     plot(tinter(:,i+index-1), datainter(:,i+index-1)./(max(abs(datainter(:,i+index-1)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
% end
% % add zeros
% % calculate the whole time window
% tstart = min(tinter(1,index:end));
% tend = max(tinter(end,index:end));
% ntotal = round((tend-tstart)/newdt)+1;
% nadd = ntotal-newnpts;
% tadd1 = zeros(nsum,1);
% nadd1 = zeros(nsum,1);
% nadd2 = zeros(nsum,1);
% for i = 1:nsum
%     tadd1(i) = tinter(1,i+index-1)-tstart;
%     nadd1(i) = round(tadd1(i)/newdt);
%     nadd2(i) = nadd-nadd1(i);
% end
% % front and behind time window added on the original time trace
% figure (3)
% tfull = zeros(ntotal,nsum);
% datafull = zeros(ntotal,nsum);
% for i=1:nsum
%     tfull1 = zeros(nadd1(i),1);
%     tfull2 = zeros(nadd2(i),1);
%     for j = 1:nadd1(i)
%         tfull1(j)=tinter(1,i+index-1)-newdt*(nadd1(i)+1-j);
%     end
%     for j = 1:nadd2(i)
%         tfull2(j)=tinter(end,i+index-1)+newdt*j;
%     end
%     % concatenate time and data arrays
%     tfull(:,i) = cat(1,tfull1,tinter(:,i+index-1),tfull2);
%     datafull(:,i) = cat(1,zeros(nadd1(i),1),datainter(:,i+index-1),zeros(nadd2(i),1));
%     plot(tfull(:,i), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
% end

%% 2. interpolation, after test, more smooth then resample
% newdt = 1/100; 
% for i = 1:sa
%     tinter(:,i) = t(1,i):newdt:t(end,i);
%     datainter(:,i) = interp1(t(:,i), data(:,i), tinter(:,i), 'spline');
% end
% % shift and plot
% for i=1:sa
%     if dist(i) > 300
%         index = i;
%         break
%     end
% end
% [newnpts, ttt] = size(tinter(:,1));       % ttt is useless
% nsum = sa-index+1;
% tshift = zeros(nsum,1);
% nshift = zeros(nsum,1);
% for i = 1:nsum
%     tshift(i) = roundn(dist(i+index-1)./4.756, -2);     % save two digits of decimal
%     nshift(i) = round(tshift(i)/newdt);
% end
% figure (2)
% for i = 1:nsum
%     tinter(:,i+index-1) =  tinter(:,i+index-1) - tshift(i);
%     plot(tinter(:,i+index-1), datainter(:,i+index-1)./(max(abs(datainter(:,i+index-1)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
% end
% % add zeros
% % calculate the whole time window
% tstart = min(tinter(1,index:end));
% tend = max(tinter(end,index:end));
% ntotal = round((tend-tstart)/newdt)+1;
% nadd = ntotal-newnpts;
% tadd1 = zeros(nsum,1);
% nadd1 = zeros(nsum,1);
% nadd2 = zeros(nsum,1);
% for i = 1:nsum
%     tadd1(i) = tinter(1,i+index-1)-tstart;
%     nadd1(i) = round(tadd1(i)/newdt);
%     nadd2(i) = nadd-nadd1(i);
% end
% % front and behind time window added on the original time trace
% figure (3)
% tfull = zeros(ntotal,nsum);
% datafull = zeros(ntotal,nsum);
% for i=1:nsum
%     tfull1 = zeros(nadd1(i),1);
%     tfull2 = zeros(nadd2(i),1);
%     for j = 1:nadd1(i)
%         tfull1(j)=tinter(1,i+index-1)-newdt*(nadd1(i)+1-j);
%     end
%     for j = 1:nadd2(i)
%         tfull2(j)=tinter(end,i+index-1)+newdt*j;
%     end
%     % concatenate time and data arrays
%     tfull(:,i) = cat(1,tfull1,tinter(:,i+index-1),tfull2);
%     datafull(:,i) = cat(1,zeros(nadd1(i),1),datainter(:,i+index-1),zeros(nadd2(i),1));
%     plot(tfull(:,i), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
% end

%% 3. directly shift in frequency domain, use specshift
for i=1:sa
    if dist(i) > 300
        index = i;
        break
    end
end
nsum = sa-index+1;
tshift = zeros(nsum,1);
nshift = zeros(nsum,1);
for i = 1:nsum
    tshift(i) = t(1,index)-t(1,i+index-1)+dist(i+index-1)./4.756;     
    nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed
    
end
datafull=[];
for i = 1:nsum
    dout = specshift(data(:,i+index-1),nshift(i)-nshift(1));     % shift according to the first trace
    datafull = [datafull dout];      % all data has the same length as original
end
for i = 1:nsum
    tfull(:,i) = t(:,index)-tshift(1);      % all data has the same time trace, i.e. first time trace - tshift
end
figure (2)
for i = 1:nsum
    plot(tfull(:,1), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end

%% interpolate maybe
newdt = 0.1;
for i=1:nsum
    tinter(:,i) = tfull(1,1): newdt: tfull(end,1);
    datainter(:,i) =  interp1(tfull(:,i), datafull(:,i), tinter(:,i), 'spline');
end
figure (3)
for i = 1:nsum
   plot(tinter(:,i), datainter(:,i)./(max(abs(datainter(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
% set(gca,'XLim',[0 400]);  % 设置x轴的间隔，长度
set(gca,'XTick',0:2:30);

%% radon transform
% tau = T - PX
% to make sure dat can be aligned more precisely, data has been
% interpoalate to a very high rate, but it will consume much time to run
% radon transform then, so may be it is a good choice to recover the
% original sampling rate.

vs = [3.04; 3.51; 3.77; 4.756];
vp = [5.09208; 5.90612; 6.35644; 8.06419];
% if shift already, preferring phase q is about 0, range should contain negative ones
% if not, q may be positive.
q = -0.2:0.001:0.2;       % dp more small, more smooth
flag = 1;     % linear transfrom
%q = -0.044: 0.001: 0.044; 
%flag = 2;        % parabola transfrom
distuse = (dist(index:end)'-300) ;   % , ==>tau=T-P(X-300), turn distance in km to m
noise = 0.5;      % big noise is better than small
% if interpolate, use datainter, tinter, newdt;   else use datafull, tinter, dt
m = f_radon(datainter(700:950,:) , newdt , distuse , q , flag , noise);
figure (4)
imagesc(q,tinter(700:950,1),m);
%colormap(hot); 
colorbar;
caxis([-5e-4,5e-4]);
% minm = min(min(m(450,)));
% [row, col]=find(m==minm);
% figure (4)
% imagesc(q,tinter(400:1400,1),m(400:1400,:));
% colorbar ;