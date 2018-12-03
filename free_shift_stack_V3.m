% free_shift_stack_V3.0
%
% align sac data acording to different ray parameter P 
% eg. Sn, sSn, Sb, Pn, pPn, Pb
% after shift, sum all to enhance the spefic phase
% Difference from v2.0:
%     change the polarity only once before use different P parameters,
%     instead of changing them in the loop
%
% Author: C. Song,  2017.4.21
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
y = (dist - 600)./4.5+145;
%plot(y, dist, 'k-', 'LineWidth',1.5);

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
mindist = 490;
for i=1:sa
    if dist(i) > mindist           % sSn appears after 300 km; Pn > 270km; pPn > 450km ; Sn > ? km
        index = i;
        break
    end
end
nsum = sa-index+1;
tshift = zeros(nsum,1);
nshift = zeros(nsum,1);

% select a reference distance
avedist = (min(dist(index:end)) + max(dist))/2;
[~, ind4] = min(abs(dist-avedist));
refdist = dist(ind4);

%% 1. use a rough velocity to shift data and change polarity of original data
% shift according to a rough velocity
refvel = 6.5;
for i = 1:nsum
    tshift(i) = t(1, ind4)-t(1, i+index-1)+dist(i+index-1)./refvel;         % 8.20 km/s
    nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
end
data1 = zeros(sc, nsum);
for i = 1:nsum
    dout = specshift(data(:, i+index-1), nshift(i)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
    data1(:, i) = dout;      % all data has the same length as original
end
t1= zeros(sc, nsum);
for i = 1:nsum
    t1(:, i) = t(:, ind4);      % all data has the same time trace, i.e. first time trace - tshift
end
figure (2)
for i = 1:nsum
    plot(t1(:,1), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:200);

% change polarity to original data

tstart = 86;      % this time window contains Pn - pPn
tend = 90.5;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;
data2 = zeros(sc, nsum);
for i = 1:nsum
        % use sign of max or min value within the tw to make negative to be
        % positive, without normalization
        [~, ind1] = max(abs(data1(pstart:pend, i)));
        data2(:, i) = data(:, i+index-1).*sign(data1(pstart+ind1-1, i));    % use data1 get polarity, correct data, get data2
end
t2 = t(:, index:end);

%% 2. use velocity first
% set velocity range
vpn = 5.0: 0.01: 8.0;        % 理论的Sn、sSn速度是4.756；Pn是8.06419  
nvpn = length(vpn);
datasum = zeros(ntw, nvpn);
tuse = zeros(ntw, nvpn);
figure (3)
for j = 1:nvpn
    j
%      j=61;
    for i = 1:nsum
        %tshift(i) = dist(i+index-1)./vpn(j);
        tshift(i) = t2(1, 1+ind4-index)-t2(1, i)+dist(i+index-1)./vpn(j);    
        nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
    end
    data3 = zeros(sc, nsum);
    for i = 1:nsum
        dout = specshift(data2(:, i), nshift(i)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
        data3(:, i) = dout;      % all data has the same length as original
    end
    t3 = zeros(sc, nsum); 
    for i = 1:nsum
        t3(:, i) = t2(:, 1+ind4-index);      % all data has the same time trace, i.e. first time trace - tshift
    end
%     figure (2)
%     for i = 1:nsum
%         plot(tfull(:,1), datafull(:,i)./(max(abs(datafull(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
%     end
%     set(gca,'XTick',0:2:100);

    % cut window for sum
    datause = data3(pstart:pend, :);    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = t3(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j)*100 + j, 'LineWidth',1);hold on;
end

figure (4)
imagesc(tuse(:, 1), vpn, datasum');
colormap(jet);
colorbar;
% caxis([-0.06, 0.06]);
xlabel('time');
ylabel('velocity');

% get vpn, taupn
dsumpn = zeros(2, nvpn);
tpns=71;
tpne=73;
[~, ppns] = min(abs((tuse(:, 1)-tpns)));
ppne = ppns+round((tpne-tpns)./dt);
for i = 1:nvpn
    [dsumpn(1, i), dsumpn(2, i)] = max(datasum(ppns:ppne, i), [], 1);
end
[~, ind2] = max(dsumpn(1,:));
bestvpn = vpn(ind2);
ind3 = dsumpn(2, ind2);
besttpn = tuse(ppns+ind3-1, 1);
besttaupn= tuse(ppns+ind3-1, 1) - dist(ind4)/bestvpn;

% get vppn, tauppn
dsumppn = zeros(2, nvpn);
tppns=76;
tppne=79;
[~, pppns] = min(abs((tuse(:, 1)-tppns)));
pppne = pppns+round((tppne-tppns)./dt);
for i = 1:nvpn
    [dsumppn(1, i), dsumppn(2, i)] = max(datasum(pppns:pppne, i), [], 1);
end
[~, ind2] = max(dsumppn(1,:));
bestvppn = vpn(ind2);
ind3 = dsumppn(2, ind2);
besttppn = tuse(ppns+ind3-1, 1);
besttauppn= tuse(pppns+ind3-1, 1) - dist(ind4)/bestvppn;
dtau = besttppn - besttpn;


%% 3. prepare data for shift by ray parameter P 
% previously we obtain data1 by shift data using 8.20km/s 

% change polarity to data1
data4 = zeros(sc, nsum);
for i = 1:nsum
        % use sign of max or min value within the tw to make negative to be
        % positive, without normalization
        [~, ind1] = max(abs(data1(pstart:pend, i)));
        data4(:, i) = data1(:, i).*sign(data1(pstart+ind1-1, i));    % use data1 get polarity, correct data, get data2
end
t4 = t1;
figure (8)
for i = 1:nsum
    %ntw = pend-pstart+1;
    temp1 = data4(:, i)./(max(abs(data4(:, i)))) .*15.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t4(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.53 0.81 0.98], 'LineWidth',1); hold on;
    fill(t4(:, i), temp2+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
%     plot(t2(figpstart:figpend, i), temp1(figpstart:figpend)+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.39 0.58 1], 'LineWidth',1); hold on;
%     fill(t2(figpstart:figpend, i), temp2(figpstart:figpend)+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
%fill([t(:, 1), fliplr(t(:, 1))], [y, fliplr(a)], 'r', 'edgealpha', 0);
end
set(gca, 'XTick', 50:10:100);
set(gca, 'xlim', [50, 100]);
set(gca, 'ylim', [470, 620]);

% suppress the data4 from Pg to end 
[~, aa] = min(abs(t4(:, 1)- tend));
data4(aa:end, :) = data4(aa:end, :)*0.001;       % suppress by set the amp to 1%  
% Only uncomment following when using t comp. data, suppress data4 from begining to Sn as well
[~, aa] = min(abs(t4(:, 1)- tstart));
data4(1:aa, :) = data4(1:aa, :)*0.001;
% ends here
figure (5)
for i = 1:nsum
    plot(t4(:, 1), data4(:,i)./(max(abs(data4(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:100);

%% 4. use ray parameter P
% set parameter range
ppn = -0.06:0.001:0.06;
nppn = length(ppn);
% time window for P need to be more wide
tstart = 82; 
tend = 95;
[~, pstart] = min(abs((t4(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;
datasum = zeros(ntw, nppn);
tuse = zeros(ntw, nppn);
figure (6)
for j = 1:nppn
     j 
    for i = 1:nsum
        tshift(i) = dist(i+index-1).*ppn(j);
        nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
    end
    data5 = zeros(sc, nsum);
    for i = 1:nsum
        dout = specshift(data4(:, i), nshift(i)-nshift(1+ind4-index));     % shift according to the first trace
        data5(:, i) = dout;      % all data has the same length as original
    end
    t5 = t4;
%     figure (2)
%     for i = 1:nsum
%         plot(tf(:,1), dataf(:,i)./(max(abs(dataf(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
%     end
%     set(gca,'XTick',-50:2:100);

    datause = data5(pstart:pend, :);    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = t5(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j)*100 + j, 'LineWidth',1);hold on;
end

figure (7)
imagesc(tuse(:, 1), ppn, datasum');
colormap(jet);
colorbar;
% caxis([-0.06, 0.06]);
xlabel('time');
ylabel('ray parameter');

% get parameter of Pn, taupn
dsumpn = zeros(2, nppn);
tpns=71;
tpne=75;
[~, ppns] = min(abs((tuse(:, 1)-tpns)));
ppne = ppns+round((tpne-tpns)./dt);
for i = 1:nppn
    [dsumpn(1, i), dsumpn(2, i)] = max(datasum(ppns:ppne, i), [], 1);
end
[~, ind2] = max(dsumpn(1,:));
bestppn = 1/8.2+ppn(ind2);
ind3 = dsumpn(2, ind2);
besttpn = tuse(ppns+ind3-1, 1);
besttaupn= tuse(ppns+ind3-1, 1) - dist(ind4)*bestppn;

% get parameter of pPn, tauppn
dsumppn = zeros(2, nppn);
tppns=76;
tppne=79;
[~, pppns] = min(abs((tuse(:, 1)-tppns)));
pppne = pppns+round((tppne-tppns)./dt);
for i = 1:nppn
    [dsumppn(1, i), dsumppn(2, i)] = max(datasum(pppns:pppne, i), [], 1);
end
[~, ind2] = max(dsumppn(1,:));
bestpppn = 1/8.2+ppn(ind2);
ind3 = dsumppn(2, ind2);
besttppn = tuse(pppns+ind3-1, 1);
besttauppn= tuse(pppns+ind3-1, 1) - dist(ind4)*bestpppn;
dtau = besttppn - besttpn;


% %regmax = imregionalmax(datasum);
% [locmax, indices] = localmax(datasum, [], false);
% [maxval(:, 1), maxval(:, 2)] = find(locmax);            %1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(maxval(:, 1));
% %maxval = zeros(dim,1);
% for i = 1:dim
%     maxval(i, 3) = tuse(maxval(i, 1));
%     maxval(i, 4) = ppn(maxval(i, 2));
% %     maxval(i, 4) = 1/refvel+ppn(maxval(i, 2));
%     maxval(i, 5) = datasum(maxval(i, 1), maxval(i, 2));         %3rd col = max. value of regional max. of data
% end
