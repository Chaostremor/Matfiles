% align sac data acording to different ray parameter P 
% eg. Sn, Sg, Sb, Pn, Pb
% after shift, sum all to enhance the spefic phase
%
% Author: C. Song,  2017.4.8
% 
clear ; clc ; close all;
%% read data, real or synthetic
%datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
%datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录
datadir = 'G:\Alxa\nodecimate\3test600\' ;    % 数据所在目录
%fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接
fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
%filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data = zeros(sc,sa);
dist = zeros(sa,1);       % dist是头段中的震中距
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end
% plot original data
figure (1)
for i = 1:sa
   % normalize
   %plot(t(:,i), data(:,i)*100 + double(dist(i)), 'LineWidth',1);hold on;'r-'
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)) , 'LineWidth',1);hold on;
   %data(:,i) = data(:,i)./(max(abs(data(:,i))));        % normalize data
end

% 
% for i = 1:sa
%     filename = strcat(datadir,strcat(stnm(i,:),'.r')) ;  % 获取文件路径
%     %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
%     [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
%     dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
% end
% for i = 1:sa
%    % normalize
%    %plot(t(:,i), data(:,i)*100 + double(dist(i)), 'LineWidth',1);hold on;
%    plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'b-', 'LineWidth',1);hold on;
%    %data(:,i) = data(:,i)./(max(abs(data(:,i))));        % normalize data
% end
% set(gca,'XTick',0:5:round(max(t(end,:))));
% set(gca,'YTick',round(min(dist))-10:10:round(max(dist))+10);


%% Here only use 3rd way to shift the data according to one specific parameter p( or velocity )
%% 3. directly shift in frequency domain, use specshift
% select distance
mindist = 400;
for i=1:sa
    if dist(i) > mindist           % sSn appears after 300 km; Pn > 270km; pPn > 450km 
        index = i;
        break
    end
end
nsum = sa-index+1;
tshift = zeros(nsum,1);
nshift = zeros(nsum,1);

%% use velocity first
% set velocity range
vpn = 7.7: 0.01: 8.3;        % 理论的Sn、sSn速度是4.756；Pn是8.06419  
nvpn = length(vpn);
%datasum = zeros(41, nvsn);
figure (3)
for j = 1:nvpn
     %j=31;
    for i = 1:nsum
        tshift(i) = t(1, index)-t(1, i+index-1)+dist(i+index-1)./vpn(j);    
        nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
    end
    datafull = zeros(sc, nsum);
    for i = 1:nsum
        dout = specshift(data(:, i+index-1), nshift(i)-nshift(1));     % i(ref)=1+ind4-index, shift according to the reference trace
        %dout = dout./max(abs(dout));        % normalize
        %[~, ind4] = max(abs(dout));
        %dout = dout./dout(ind4);
        datafull(:, i) = dout;      % all data has the same length as original
    end
    tfull = zeros(sc, nsum); 
    for i = 1:nsum
        tfull(:, i) = t(:, index)-tshift(1);      % all data has the same time trace, i.e. first time trace - tshift
    end
%     figure (2)
%     for i = 1:nsum
%         plot(tfull(:,1), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
%     end
%     set(gca,'XTick',0:2:100);
    % cut window for normalization and sum
%     tstart = 6; 
%     tend = 18.5;
    tstart = 6 + 0.5*(dist(index)+dist(end))./vpn(1) - 0.5*(dist(index)+dist(end))./vpn(j);
    tend = 15+ 0.5*(dist(index)+dist(end))./vpn(1) - 0.5*(dist(index)+dist(end))./vpn(j);
    %pstart = find(abs(tfull(:, 1)-tstart)==min(abs((tfull(:,1)-tstart))));   % two ways to get start point
    [~, pstart] = min(abs((tfull(:,1)-tstart)));
    pend = pstart+round((tend-tstart)./dt);
    ntw = pend-pstart+1;
    datause = zeros(ntw, nsum);
    for i = 1:nsum
        %datause(:, i) = datafull(pstart:pend, i);     % no normalize
        %datause(:, i) = datafull(pstart:pend, i)./(max(abs(datafull(pstart:pend, i))));   % normalize
        % use max or min value within the tw to normalize, AND make negative to be positive
        [~, ind1] = max(abs(datafull(pstart:pend, i)));
        datause(:, i) = datafull(pstart:pend, i)./datafull(pstart+ind1-1, i);
        % use sign of max or min value within the tw to make negative to be positive
%         [~, ind1] = max(abs(datafull(pstart:pend, i)));
%         datause(:, i) = datafull(pstart:pend, i).*sign(datafull(pstart+ind1-1, i));
    end
    %datasum(:, j) = sum(abs(datause),2);       % direct sum absolute value    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = tfull(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j) + j, 'LineWidth',1);hold on;
%     plot(tfull(pstart:pend, 1), datasum(:, j) + j, 'LineWidth',1);hold on;
%     plot(tfull(pstart:pend, 1), datasum(:, j)*100 + j, 'LineWidth',1);hold on;   
end

% add zeros, to make the seq. same length 
% calculate the whole time window
tfstart = min(tuse(1, :));      
tfend = max(tuse(end, :));
ntotal = round((tfend-tfstart)/dt)+1;
nadd = ntotal-ntw;
tadd1 = zeros(nvpn,1);     % differential time before current seq.
nadd1 = zeros(nvpn,1);    % num. of points before
nadd2 = zeros(nvpn,1);    % num. of points after     
for i = 1:nvpn
    tadd1(i) = tuse(1, i)-tfstart;
    nadd1(i) = round(tadd1(i)/dt);
    nadd2(i) = nadd-nadd1(i);
end
% front and behind time window added on the original time trace
figure (4)
tf = zeros(ntotal, nvpn);           % full time seq.
dataf = zeros(ntotal, nvpn);     % full data seq.
for i=1:nvpn
    tf1 = zeros(nadd1(i),1);        % front segment of t seq.
    tf2 = zeros(nadd2(i),1);        % behind segment of t seq.
    for j = 1:nadd1(i)                     
        tf1(j)=tuse(1, i)-dt*(nadd1(i)+1-j);    
    end
    for j = 1:nadd2(i)     
        tf2(j)=tuse(end,i)+dt*j;
    end
    % concatenate time and data arrays
    tf(:, i) = cat(1, tf1, tuse(:, i), tf2);
    dataf(:, i) = cat(1, zeros(nadd1(i), 1), datasum(:, i), zeros(nadd2(i),1));
    plot(tf(:,i), dataf(:,i) + i, 'LineWidth',1);hold on;
end

figure (5)
imagesc(tf(:, 1), vpn, dataf');
colormap(jet);
colorbar;
% caxis([10,30]);
xlabel('tau');
ylabel('velocity');

% figure (5)
% for i = 1:nvpn
%     imagesc(tuse(:, i), vpn(i), datasum(:, i)'); hold on;
% end
% colormap(jet);
% colorbar;
% caxis([5,30]);
dsummax = zeros(2, nvpn);
for i = 1:nvpn
    [dsummax(1, i), dsummax(2, i)] = max(datasum(:, i), [], 1);
end
[bestsum, ind2] = max(dsummax(1,:));
bestvpn = vpn(ind2);
tstart = 6 + 0.5*(dist(index)+dist(end))./vpn(1) - 0.5*(dist(index)+dist(end))./vpn(ind2);
tend = 10.5 + 0.5*(dist(index)+dist(end))./vpn(1) - 0.5*(dist(index)+dist(end))./vpn(ind2);
[~, pstart] = min(abs((tfull(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ind3 = dsummax(2, ind2);
besttau = tfull(pstart+ind3-1, 1);

%% get the result: best vpn be around 8 km/s, then shift the whole data according to this velocity
for i = 1:nsum
    tshift(i) = t(1, index)-t(1, i+index-1)+dist(i+index-1)./8.0;     % use the 1st trace >300 as the reference
    nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
end
datap = zeros(sc, nsum);
for i = 1:nsum
    dout = specshift(data(:,i+index-1), nshift(i)-nshift(1));     % shift according to the first trace
    datap(:, i) = dout;      % all data has the same length as original
end
tp = zeros(sc, nsum);
for i = 1:nsum
    tp(:, i) = t(:, index)-tshift(1);      % all data has the same time trace, i.e. first time trace - tshift
end
figure (6)
for i = 1:nsum
    plot(tp(:,1), datap(:,i)./(max(abs(datap(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:50);

%% use ray parameter P
% set parameter range
ppn = -0.02:0.001:0.02;
nppn = length(ppn);
datasum = [];
tuse = [];
figure (7)
for j = 1:nppn
%     j = 41;
    for i = 1:nsum
        tshift(i) = dist(i+index-1).*ppn(j);
        nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
    end
    datafull = zeros(sc, nsum);
    for i = 1:nsum
        dout = specshift(datap(:, i), nshift(i)-nshift(1));     % shift according to the first trace
        datafull(:, i) = dout;      % all data has the same length as original
    end
    tfull = zeros(sc, nsum); 
    for i = 1:nsum
        tfull(:, i) = tp(:, 1)-tshift(1);      % all data has the same time trace, i.e. first time trace - tshift
    end
%     figure (2)
%     for i = 1:nsum
%         plot(tfull(:,1), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
%     end
%     set(gca,'XTick',-50:2:100);
    % cut window for normalization and sum
    tstart = 0; 
    tend = 40;
%     tstart = 22 + 0.5*(dist(index)+dist(end)).*( 0.05- ppn(j));
%     tend = 27 + 0.5*(dist(index)+dist(end)).*(ppn(1) - ppn(j));
    %pstart = find(abs(tfull(:, 1)-tstart)==min(abs((tfull(:,1)-tstart))));   % two ways to get start point
    [~, pstart] = min(abs((tfull(:,1)-tstart)));
    pend = pstart+round((tend-tstart)./dt);
    ntw = pend-pstart+1;
    datause = zeros(ntw, nsum);
    for i = 1:nsum
        %datause(:, i) = datafull(pstart:pend, i);     % no normalize
        % use max or min value within the tw to normalize, AND make negative to be positive
%         [~, ind1] = max(abs(datafull(pstart:pend, i)));
%         datause(:, i) = datafull(pstart:pend, i)./datafull(pstart+ind1-1, i);
        % use sign of max or min value within the tw to make negative to be positive
        [~, ind1] = max(abs(datafull(pstart:pend, i)));
        datause(:, i) = datafull(pstart:pend, i).*sign(datafull(pstart+ind1-1, i));
    end
    %datasum(:, j) = sum(abs(datause),2);       % direct sum absolute value    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = tfull(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j)*10 + j, 'LineWidth',1);hold on;
%     plot(tfull(pstart:pend, 1), datasum(:, j) + j, 'LineWidth',1);hold on;
end

% % add zeros, to make the seq. same length 
% % calculate the whole time window
% tfstart = min(tuse(1, :));      
% tfend = max(tuse(end, :));
% ntotal = round((tfend-tfstart)/dt)+1;
% nadd = ntotal-ntw;
% tadd1 = zeros(nvpn,1);     % differential time before current seq.
% nadd1 = zeros(nvpn,1);    % num. of points before
% nadd2 = zeros(nvpn,1);    % num. of points after     
% for i = 1:nvpn
%     tadd1(i) = tuse(1, i)-tfstart;
%     nadd1(i) = round(tadd1(i)/dt);
%     nadd2(i) = nadd-nadd1(i);
% end
% % front and behind time window added on the original time trace
% figure (4)
% tf = zeros(ntotal, nvpn);           % full time seq.
% dataf = zeros(ntotal, nvpn);     % full data seq.
% for i=1:nvpn
%     tf1 = zeros(nadd1(i),1);        % front segment of t seq.
%     tf2 = zeros(nadd2(i),1);        % behind segment of t seq.
%     for j = 1:nadd1(i)                     
%         tf1(j)=tuse(1, i)-dt*(nadd1(i)+1-j);    
%     end
%     for j = 1:nadd2(i)     
%         tf2(j)=tuse(end,i)+dt*j;
%     end
%     % concatenate time and data arrays
%     tf(:, i) = cat(1, tf1, tuse(:, i), tf2);
%     dataf(:, i) = cat(1, zeros(nadd1(i), 1), datasum(:, i), zeros(nadd2(i),1));
%     plot(tf(:,i), dataf(:,i) + i, 'LineWidth',1);hold on;
% end

figure (8)
imagesc(tuse(:, 1), ppn, datasum');
colormap(jet);
colorbar;
% caxis([10,30]);
xlabel('tau');
ylabel('ray parameter');
