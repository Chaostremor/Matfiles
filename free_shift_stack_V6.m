% free shift stack V6
% 
% the simpliest version 
%
% USAGE:
%     1. read data
%     2. unify the polarity of data
%     3. use different vel. to phrase weighed stack (PWS) all data
%     4. use different ray parameter to PWS all data
%
% Difference from V5:
%      1. no shift before ray parameter stack
%      2. no suppress to data
%      3. no normalization to sum
%      ...
%
%  Author: C.Song, 2017.5.3

clear ; clc ; close all;
%% 1. read data
datadir = 'G:\Alxa\nodecimate\3test400\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
ii=1;
filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data = zeros(sc,sa);
dist = zeros(sa,1);       % dist是头段中的震中距
pg = zeros(sa,1);
sg = zeros(sa,1);
for ii = 1:sa
    filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % 获取文件路径
    [t(:,ii), data(:,ii), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(ii) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    pg(ii) = SAChdr.times.t2;
    sg(ii) = SAChdr.times.t3;
end
% normalize to plot original data
figure
for ii = 1:sa
   plot(t(:,ii), data(:,ii)./(max(abs(data(:,ii)))) .*15.0 + double(dist(ii)) , 'LineWidth',1);hold on;
end
% y = (dist - 600)./4.5+145;
% plot(y, dist, 'k-', 'LineWidth',1.5);


%% 2. use a rough velocity to shift data and change polarity of original data
% select distance
mindist = 250;
for ii=1:sa
    if dist(ii) > mindist         % sSn appears after 300 km; Pn > 270km; pPn > 450km ; Sn > ? km
        index = ii;
        break
    end
end
nsum = sa-index+1;
tshiftp = zeros(nsum,1);
nshiftp = zeros(nsum,1);

% select a reference distance
avedist = (min(dist(index:end)) + max(dist))/2;
[~, ind4] = min(abs(dist-avedist));
refdist = dist(ind4);

% shift data and the phase of it according to a rough velocity
refvel = 8.2;                          % the reference velocity is almost Pn-pPn,  or Sn-sSn
for ii = 1:nsum
    tshiftp(ii) = t(1, ind4)-t(1, ii+index-1)+dist(ii+index-1)./refvel;         % 8.20 km/s
    nshiftp(ii) =tshiftp(ii)/dt;      % in freq. domain, non-integer is allowed    
end
data1 = zeros(sc, nsum);
for ii = 1:nsum
    dout = specshift(data(:, ii+index-1), nshiftp(ii)-nshiftp(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
    data1(:, ii) = dout;      % all data has the same length as original
end
t1= zeros(sc, nsum);
for ii = 1:nsum
    t1(:, ii) = t(:, ind4);      % all data has the same time trace, i.e. first time trace - tshift
end
figure 
for ii = 1:nsum
    plot(t1(:,1), data1(:,ii)./(max(abs(data1(:,ii)))) .*15.0 + double(dist(ii+index-1)), 'LineWidth',1);hold on;
end
plot(50, refdist, 'k.', 'MarkerSize', 10); 
set(gca,'XTick',0:2:200);

% change polarity to original data and shifted data 
tstart = 48;      % this time window ONLY contains Pn-pPn, CHANGE when ref. dist change 
tend = 52;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
data2 = zeros(sc, nsum);
phase2 = zeros(sc, nsum);
pn = zeros(nsum, 1);

% % if deal with SH comp. , still use the polarity from P-SV comp. 
% START : uncomment when dealing with SH comp.
for ii = 1:sa
    filename = strcat(datadir,strcat(stnm(ii,:),'.t')) ;  % 获取文件路径
    [~, data(:,ii), ~] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
end
% normalize to plot original data
figure
for ii = 1:sa
   plot(t(:,ii), data(:,ii)./(max(abs(data(:,ii)))) .*15.0 + double(dist(ii)) , 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:200);

refvel2 = 4.4;
tshifts =  zeros(nsum,1);
nshifts = zeros(nsum,1);
sn = zeros(nsum, 1);
for ii = 1:nsum
    tshifts(ii) = t(1, ind4)-t(1, ii+index-1)+dist(ii+index-1)./refvel2;        
    nshifts(ii) =tshifts(ii)/dt;      % in freq. domain, non-integer is allowed    
end
dataxx = zeros(sc, nsum);
for ii = 1:nsum
    dout = specshift(data(:, ii+index-1), nshifts(ii)-nshifts(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
    dataxx(:, ii) = dout;      % all data has the same length as original
end
txx= zeros(sc, nsum);
for ii = 1:nsum
    txx(:, ii) = t(:, ind4);      % all data has the same time trace, i.e. first time trace - tshift
end
figure 
for ii = 1:nsum    
    plot(txx(:,1), dataxx(:,ii)./(max(abs(dataxx(:,ii)))) .*15.0 + double(dist(ii+index-1)), 'LineWidth',1);hold on;
end
plot(50, refdist, 'k.', 'MarkerSize', 10); 
set(gca,'XTick',0:2:200);
tstarts = 81;      % this time window ONLY contains Pn-pPn, CHANGE when ref. dist change 
tends = 88;
[~, pstarts] = min(abs((txx(:,1)-tstarts)));
pends = pstarts+round((tends-tstarts)./dt);
for ii = 1:nsum    
    [~, inds] = max(abs(dataxx(pstarts:pends, ii)));        
    sn(ii) = txx(pstarts+inds-1, ii) + tshifts(ii)-tshifts(1+ind4-index);
end
% END

for ii = 1:nsum
        % use sign of max or min value within the tw to make negative to be
        % positive, without normalization
        [~, ind1] = max(abs(data1(pstart:pend, ii)));        
        pn(ii) = t1(pstart+ind1-1, ii) + tshiftp(ii)-tshiftp(1+ind4-index);
        data2(:, ii) = data(:, ii+index-1).*sign(data1(pstart+ind1-1, ii));    % use data1 get polarity, correct used data, get data2
        data1(:, ii) = data1(:, ii).*sign(data1(pstart+ind1-1, ii));      % correct data1's own polarity
        phase2(:, ii) = phase(hilbert(data2(:, ii)));                % do hilbert transform to data, and get the phase phi of the complex DATA
end
t2 = t(:, index:end);
figure
for ii = 1:nsum
    temp1 = data2(:, ii)./(max(abs(data2(:, ii)))) .*15.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t2(:, ii), temp1+ double(dist(ii+index-1)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
    fill(t2(:, ii), temp2+ double(dist(ii+index-1)), 'r', 'edgealpha', 0); hold on;
end
set(gca, 'XTick', 0:5:150);
set(gca, 'xlim', [0, 150]);
set(gca, 'ylim', [0, 420]);

% save('pgsgpicktime.mat', 'pg', 'sg', 'dist', 'pn', 'sn', 'index');

%% 3. use velocity first
% set velocity range
vel = 2.5: 0.01: 6;        % 理论的Sn、sSn速度是4.756；Pn是8.06419  
nvel = length(vel);

% cut window for sum
tstart = 76;      % in fact, it is tw contains phases at the reference distance !!!
tend = 98;       % CHANGE if the ref. dist changes
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw1 = pend-pstart+1;

% shift and sum
datasum1 = zeros(ntw1, nvel);
phasesum1 = zeros(ntw1, nvel);
pwssum1 = zeros(ntw1, nvel);
tuse1 = zeros(ntw1, nvel);
tshift =  zeros(nsum,1);
nshift = zeros(nsum,1);
power = 1;
%figure
for jj = 1:nvel
    jj
%      jj=61;
    for ii = 1:nsum
        tshift(ii) = t2(1, 1+ind4-index)-t2(1, ii)+dist(ii+index-1)./vel(jj);    
        nshift(ii) =tshift(ii)/dt;      % in freq. domain, non-integer is allowed    
    end
    data3 = zeros(sc, nsum);
    phase3 = zeros(sc, nsum);
    for ii = 1:nsum
        dout = specshift(data2(:, ii), nshift(ii)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
        data3(:, ii) = dout;      % all data has the same length as original
        pout = specshift(phase2(:, ii), nshift(ii)-nshift(1+ind4-index));
        phase3(:, ii) = pout;
    end
    t3 = zeros(sc, nsum); 
    for ii = 1:nsum
        t3(:, ii) = t2(:, 1+ind4-index);      % all data has the same time trace, i.e. first time trace - tshift
    end
    
    datause = data3(pstart:pend, :);
    datasum1(:, jj) = sum(datause, 2)./nsum;      % average linear stacked amplitude
    phaseuse = phase3(pstart:pend, :);
    phasesum1(:, jj) = abs(sum(exp(1i.*phaseuse), 2)./nsum);      % average phase stack
    pwssum1(:, jj) = datasum1(:, jj).*(phasesum1(:, jj).^power);    % average phase weighed stacked amplitude
    tuse1(:, jj) = t3(pstart:pend, 1);
    %plot(tuse1(:, jj), datasum1(:, jj)*100 + jj, 'LineWidth',1);hold on;
end

% figure
% imagesc(tuse1(:, 1), vel, datasum1');
% colormap(jet);
% colorbar;
% xlabel('time');
% ylabel('velocity');

figure
imagesc(tuse1(:, 1), vel, pwssum1');
colormap(jet);
colorbar;
xlabel('time');
ylabel('velocity');

ivel=find(vel==4.48);
figure
plot(tuse1(:, ivel), pwssum1(:, ivel), 'LineWidth',1);

% regmax = imregionalmax(pwssum1);
% gmax = max(pwssum1(:));
% [irow, icol] = find(regmax);            % 1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(irow);
% %maxval = zeros(dim,1);
% maxset = [];
% for i = 1:dim
%     maxval = pwssum1(irow(i), icol(i));     % 3rd col = max. value of regional max. of data
%     if (maxval >= 0.6*gmax)
%         maxt = tuse1(irow(i));
%         maxp = vel(icol(i));
%         maxset = [maxset; irow(i) icol(i) maxt maxp maxval];
%     end
% end
% figure
% imagesc(tuse1(:, 1), vel, pwssum1'); hold on
% colormap(jet);
% colorbar;
% xlabel('time');
% ylabel('velocity');
% plot(maxset(:, 3), maxset(:, 4), 'k.', 'MarkerSize', 8); hold on
% % rectx = [85, 85, 90, 90, 85];
% % recty = [0.017, -0.036, -0.036, 0.017, 0.017];
% % plot(rectx, recty, 'linewidth', 1.5, 'color', 'w'); hold on
% [~, indice] = max(pwssum1(:));
% [row, col] = ind2sub([ntw1, nvel], indice);
% plot(tuse1(row), vel(col), 'k*', 'MarkerSize', 10); 
% [~, indice] = min(pwssum1(:));
% [row, col] = ind2sub([ntw1, nvel], indice);
% plot(tuse1(row), vel(col), 'kx', 'MarkerSize', 10); 


%% 4. use ray parameter P
% set parameter range
rayp = 0.167: 0.001: 0.400;
nrayp = length(rayp);
% time window for P parameter need to be more wide
tstart = 76;         % the tw is now the same as it in step 3 
tend = 98;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw2 = pend-pstart+1;
datause = [];
phaseuse = [];
datasum2 = zeros(ntw2, nrayp);
phasesum2 = zeros(ntw2, nrayp);
pwssum2 = zeros(ntw2, nrayp);
tuse2 = zeros(ntw2, nrayp);
power = 1;
%figure
for jj = 1:nrayp
     jj 
    for ii = 1:nsum
        tshift(ii) = t2(1, 1+ind4-index)-t2(1, ii)+dist(ii+index-1).*rayp(jj);
        nshift(ii) =tshift(ii)/dt;      % in freq. domain, non-integer is allowed    
    end
    data4 = zeros(sc, nsum);
    phase4 = zeros(sc, nsum);
    for ii = 1:nsum
        dout = specshift(data2(:, ii), nshift(ii)-nshift(1+ind4-index));     % shift according to the first trace
        data4(:, ii) = dout;      % all data has the same length as original
        pout = specshift(phase2(:, ii), nshift(ii)-nshift(1+ind4-index));
        phase4(:, ii) = pout;
    end
    t4 = zeros(sc, nsum); 
    for ii = 1:nsum
        t4(:, ii) = t2(:, 1+ind4-index);      % all data has the same time trace, i.e. first time trace - tshift
    end
    
    datause = data4(pstart:pend, :);
    datasum2(:, jj) = sum(datause, 2)./nsum;      % average linear stacked amplitude
    phaseuse = phase4(pstart:pend, :);
    phasesum2(:, jj) = abs(sum(exp(1i.*phaseuse), 2)./nsum);        % average phase stack
    pwssum2(:, jj) = datasum2(:, jj).*(phasesum2(:, jj).^power);     % average phase weighed stacked amplitude
    tuse2(:, jj) = t4(pstart:pend, 1);
    %plot(tuse2(:, jj), datasum2(:, jj)*100 + jj, 'LineWidth',1);hold on;
end

% figure
% imagesc(tuse2(:, 1), rayp, datasum2');
% colormap(jet);
% colorbar;
% xlabel('time');
% ylabel('ray parameter');

figure
imagesc(tuse2(:, 1), rayp, pwssum2');
colormap(jet);
colorbar;
xlabel('time');
ylabel('ray parameter');

% regmax = imregionalmax(pwssum2);
% gmax = max(pwssum2(:));
% [irow, icol] = find(regmax);            % 1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(irow);
% %maxval = zeros(dim,1);
% maxset = [];
% for i = 1:dim
%     maxval = pwssum2(irow(i), icol(i));     % 3rd col = max. value of regional max. of data
%     if (maxval >= 0.6*gmax)
%         maxt = tuse2(irow(i));
%         maxp = rayp(icol(i));
%         maxset = [maxset; irow(i) icol(i) maxt maxp maxval];
%     end
%            
% end
% figure
% imagesc(tuse2(:, 1), rayp, pwssum2'); hold on
% colormap(jet);
% colorbar;
% xlabel('time');
% ylabel('ray parameter');
% plot(maxset(:, 3), maxset(:, 4), 'k.', 'MarkerSize', 8); hold on
% % rectx = [85, 85, 90, 90, 85];
% % recty = [0.017, -0.036, -0.036, 0.017, 0.017];
% % plot(rectx, recty, 'linewidth', 1.5, 'color', 'w'); hold on
% [~, indice] = max(pwssum2(:));
% [row, col] = ind2sub([ntw2, nrayp], indice);
% plot(tuse2(row), rayp(col), 'k*', 'MarkerSize', 10); 
% [~, indice] = min(pwssum2(:));
% [row, col] = ind2sub([ntw2, nrayp], indice);
% plot(tuse2(row), rayp(col), 'kx', 'MarkerSize', 10); 
