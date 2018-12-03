% free_shift_stack_V4.0
%
% align sac data acording to different ray parameter P, eg. Sn, sSn, Sb, Pn, pPn, Pb
% after shift, sum all to enhance the spefic phase
%
% Difference from v3.0:
%     change the polarity only using the Pn-pPn or Sn-sSn window, and remain
%     same when stack Pn-pPn and Pb, in case one station's polarity would
%     change when stacking different phases
%
% Author: C. Song,  2017.4.27
%

clear ; clc ; close all;
%% 0. read data
datadir = 'G:\Alxa\nodecimate\3test600\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
i=1;
filename = strcat(datadir,strcat(stnm(i,:),'.z')) ;  % 获取文件路径
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
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end
% normalize to plot original data
figure
for i = 1:sa
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)) , 'LineWidth',1);hold on;
end
% y = (dist - 600)./4.5+145;
% plot(y, dist, 'k-', 'LineWidth',1.5);

%% 1. use a rough velocity to shift data and change polarity of original data
% select distance
mindist = 400;
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

% shift according to a rough velocity
refvel = 8.2;                          % the reference velocity is almost Pn-pPn,  or Sn-sSn
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
figure 
for i = 1:nsum
    plot(t1(:,1), data1(:,i)./(max(abs(data1(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:200);

% change polarity to original data and shifted data 
tstart = 68;      % this time window ONLY contains Pn-pPn or Sn-sSn
tend = 80;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
data2 = zeros(sc, nsum);
for i = 1:nsum
        % use sign of max or min value within the tw to make negative to be
        % positive, without normalization
        [~, ind1] = max(abs(data1(pstart:pend, i)));
        data2(:, i) = data(:, i+index-1).*sign(data1(pstart+ind1-1, i));    % use data1 get polarity, correct used data, get data2
        data1(:, i) = data1(:, i).*sign(data1(pstart+ind1-1, i));      % correct data1's own polarity
end
t2 = t(:, index:end);

%% 2. use velocity first
% set velocity range
vel = 6.0: 0.01: 9.0;        % 理论的Sn、sSn速度是4.756；Pn是8.06419  
nvel = length(vel);

% cut window for sum
tstart = 70;      % this time window contains the phase wanted to stacked out
tend = 85;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;

% shift and sum
datasum = zeros(ntw, nvel);
tuse = zeros(ntw, nvel);
figure
for j = 1:nvel
    j
%      j=61;
    for i = 1:nsum
        tshift(i) = t2(1, 1+ind4-index)-t2(1, i)+dist(i+index-1)./vel(j);    
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
    
    datause = data3(pstart:pend, :);    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = t3(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j)*100 + j, 'LineWidth',1);hold on;
end

figure
imagesc(tuse(:, 1), vel, datasum');
colormap(jet);
colorbar;
% caxis([-0.06, 0.06]);
xlabel('time');
ylabel('velocity');

%% 3. prepare data for shift by ray parameter P 
% NOTE: we already obtained data2 out of data, whose polarity was corrected, 
%            so we could just shift data2 according to specific phase  

% shift data2 by approximate vel. of specific phase
refvel =  6.5;                          % the reference velocity is of any specific phase 
for i = 1:nsum
    tshift(i) = t2(1, 1+ind4-index)-t2(1, i)+dist(i+index-1)./refvel;         % 8.20 km/s
    nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
end
data4 = zeros(sc, nsum);
for i = 1:nsum
    dout = specshift(data2(:, i), nshift(i)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
    data4(:, i) = dout;      % all data has the same length as original
end
t4= zeros(sc, nsum);
for i = 1:nsum
    t4(:, i) = t(:, ind4);      % all data has the same time trace, i.e. first time trace - tshift
end

figure
for i = 1:nsum
    temp1 = data4(:, i)./(max(abs(data4(:, i)))) .*15.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t4(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
    fill(t4(:, i), temp2+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
%     plot(t2(figpstart:figpend, i), temp1(figpstart:figpend)+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.39 0.58 1], 'LineWidth',1); hold on;
%     fill(t2(figpstart:figpend, i), temp2(figpstart:figpend)+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
%fill([t(:, 1), fliplr(t(:, 1))], [y, fliplr(a)], 'r', 'edgealpha', 0);
end
set(gca, 'XTick', 0:5:200);
set(gca, 'xlim', [0, 200]);
set(gca, 'ylim', [390, 620]);

% suppress the data4 from Pg/Sg to end
suptend = 84;
[~, aa] = min(abs(t4(:, 1)- suptend));
data4(aa:end, :) = data4(aa:end, :)*0.001;       % suppress by set the amp to 1%  
% % Only uncomment following when using t comp. data, suppress data4 from begining to Sn as well
% suptstart = 86;
% [~, aa] = min(abs(t4(:, 1)- suptstart));
% data4(1:aa, :) = data4(1:aa, :)*0.001;
% % ends here
figure
for i = 1:nsum
    plot(t4(:, 1), data4(:,i)./(max(abs(data4(:,i)))) .*15.0 + double(dist(i+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:100);

%% 4. use ray parameter P
% set parameter range
rayp = -0.06:0.001:0.06;
nrayp = length(rayp);
% time window for P parameter need to be more wide
tstart = 60; 
tend = 90;
[~, pstart] = min(abs((t4(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;
datasum = zeros(ntw, nrayp);
tuse = zeros(ntw, nrayp);
figure
for j = 1:nrayp
     j 
    for i = 1:nsum
        tshift(i) = dist(i+index-1).*rayp(j);
        nshift(i) =tshift(i)/dt;      % in freq. domain, non-integer is allowed    
    end
    data5 = zeros(sc, nsum);
    for i = 1:nsum
        dout = specshift(data4(:, i), nshift(i)-nshift(1+ind4-index));     % shift according to the first trace
        data5(:, i) = dout;      % all data has the same length as original
    end
    t5 = t4;

    datause = data5(pstart:pend, :);    
    datasum(:, j) = sum(datause, 2);      % direct sum
    tuse(:, j) = t5(pstart:pend, 1);
    plot(tuse(:, j), datasum(:, j)*100 + j, 'LineWidth',1);hold on;
end

% plot datasum
figure
imagesc(tuse(:, 1), rayp, datasum');
colormap(jet);
colorbar;
% caxis([-0.06, 0.06]);
xlabel('time');
ylabel('ray parameter');

% %% 5. find local maximum locations and mark them
% regmax = imregionalmax(datasum);
% [irow, icol] = find(regmax);            % 1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(irow);
% %maxval = zeros(dim,1);
% maxset = [];
% for i = 1:dim
%     maxval = datasum(irow(i), icol(i));     % 3rd col = max. value of regional max. of data
%     if (maxval >= 3e-3)
%         maxt = tuse(irow(i));
%         maxp = ppn(icol(i));
%         maxset = [maxset; irow(i) icol(i) maxt maxp maxval];
%     end
%            
% end
% figure
% imagesc(tuse(:, 1), rayp, datasum'); hold on
% colormap(jet);
% colorbar;
% % caxis([-0.06, 0.06]);
% xlabel('time');
% ylabel('ray parameter');
% plot(maxset(:, 3), maxset(:, 4), 'k.', 'MarkerSize', 8); hold on
% rectx = [85, 85, 90, 90, 85];
% recty = [0.017, -0.036, -0.036, 0.017, 0.017];
% plot(rectx, recty, 'linewidth', 1.5, 'color', 'w'); hold on
% [~, indice] = max(datasum(:));
% [row, col] = ind2sub([ntw, nppn], indice);
% plot(tuse(row), ppn(col), 'k*', 'MarkerSize', 10); 
% [~, indice] = min(datasum(:));
% [row, col] = ind2sub([ntw, nppn], indice);
% plot(tuse(row), ppn(col), 'kx', 'MarkerSize', 10); 
% 
% %% 6. plot the potential phase curves on polarity-corrected original data to check the validity
% figure
% for i = 1:nsum
%     temp1 = data2(:, i)./(max(abs(data2(:, i)))) .*15.0 ;
%     temp2 = (0.5*(temp1+abs(temp1)));
%     plot(t2(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
%     fill(t2(:, i), temp2+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
% %     plot(t2(figpstart:figpend, i), temp1(figpstart:figpend)+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.39 0.58 1], 'LineWidth',1); hold on;
% %     fill(t2(figpstart:figpend, i), temp2(figpstart:figpend)+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
% %fill([t(:, 1), fliplr(t(:, 1))], [y, fliplr(a)], 'r', 'edgealpha', 0);
% end
% set(gca, 'XTick', 20:10:200);
% set(gca, 'xlim', [20, 200]);
% set(gca, 'ylim', [470, 620]);
% 
% distcurve = 490: 0.1: 600;
% p1 = 1/7.0-0.024;         % pPn
% tau1 = 82.8;
% y1 = p1.*(distcurve - refdist )+tau1;
% plot(y1, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p2 = 0.118;                   % Pn
% tau2 = 72.3;
% y2 = p2.*(distcurve- 499.54)+tau2;     
% plot(y2, distcurve, 'linestyle', '--', 'color', 'm', 'LineWidth', 2); hold on;
% 
% p3 = 1/6.5+0.007;          % Pb
% tau3 = 89.37;
% y3 = p3.*(distcurve - refdist)+tau3;
% plot( y3, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% a1 = 5.969;
% b1 = 35.05;
% vpg = b1/a1;
% pgcurve = a1.*sqrt((distcurve./b1).^2+1);
% plot(pgcurve, distcurve, 'b-', 'LineWidth', 2); hold on;

% distcurve = 470: 0.1: 610;
% p1 = 1/5.0+0.016;         % Sn
% tau1 = 132.3;
% y1 = p1.*(distcurve - 524.2231 )+tau1;
% plot(y1, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;

% p2 = 1/5.0+0.014;                   % sSn
% tau2 = 139.8;
% y2 = p2.*(distcurve- 524.2231)+tau2;     
% plot(y2, distcurve, 'linestyle', '--', 'color', 'y', 'LineWidth', 2); hold on;
% 
% p3 = 1/4.5+0.005;          % Sb
% tau3 = 152.1;
% y3 = p3.*(distcurve - refdist)+tau3;
% plot( y3, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% a2 = 10.34;
% b2 = 35.95;
% vsg = b2/a2;
% sgcurve = a2.*sqrt((distcurve./b2).^2+1);
% plot(sgcurve, distcurve, 'b-', 'LineWidth', 2);
% 
% text(150, 610, 'Sn');
% text(156, 610, 'sSn');
% text(165, 610, 'Sb');
% text(175, 610, 'Sg');