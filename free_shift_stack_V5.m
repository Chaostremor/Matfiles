% free shift stack V5
%
% USAGE:
%     align sac data acording to different ray parameter P, eg. Sn, sSn, Sb, Pn, pPn, Pb
%     after shift, sum all to enhance the spefic phase
%
% Difference from v4.0:
%     use the phrase weighed stack method (PWS) by Schimmel and Paulssen (1997)
%
% Author:  C.Song, 2017.4.31


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
ii=1;
filename = strcat(datadir,strcat(stnm(ii,:),'.t')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data = zeros(sc,sa);
dist = zeros(sa,1);       % dist是头段中的震中距
for ii = 1:sa
    filename = strcat(datadir,strcat(stnm(ii,:),'.t')) ;  % 获取文件路径
    [t(:,ii), data(:,ii), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(ii) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end
% normalize to plot original data
figure
for ii = 1:sa
   plot(t(:,ii), data(:,ii)./(max(abs(data(:,ii)))) .*15.0 + double(dist(ii)) , 'LineWidth',1);hold on;
end
% y = (dist - 600)./4.5+145;
% plot(y, dist, 'k-', 'LineWidth',1.5);


%% 1. use a rough velocity to shift data and change polarity of original data
% select distance
mindist = 400;
for ii=1:sa
    if dist(ii) > mindist         % sSn appears after 300 km; Pn > 270km; pPn > 450km ; Sn > ? km
        index = ii;
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

% shift data and the phase of it according to a rough velocity
refvel = 5.0;                          % the reference velocity is almost Pn-pPn,  or Sn-sSn
for ii = 1:nsum
    tshift(ii) = t(1, ind4)-t(1, ii+index-1)+dist(ii+index-1)./refvel;         % 8.20 km/s
    nshift(ii) =tshift(ii)/dt;      % in freq. domain, non-integer is allowed    
end
data1 = zeros(sc, nsum);
for ii = 1:nsum
    dout = specshift(data(:, ii+index-1), nshift(ii)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
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
set(gca,'XTick',0:2:200);

% change polarity to original data and shifted data 
tstart = 120;      % this time window ONLY contains Pn-pPn or Sn-sSn
tend = 134;
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
data2 = zeros(sc, nsum);
phase2 = zeros(sc, nsum);
for ii = 1:nsum
        % use sign of max or min value within the tw to make negative to be
        % positive, without normalization
        [~, ind1] = max(abs(data1(pstart:pend, ii)));
        data2(:, ii) = data(:, ii+index-1).*sign(data1(pstart+ind1-1, ii));    % use data1 get polarity, correct used data, get data2
        data1(:, ii) = data1(:, ii).*sign(data1(pstart+ind1-1, ii));      % correct data1's own polarity
        phase2(:, ii) = phase(hilbert(data2(:, ii)));                % do hilbert transform to data, and get the phase phi of the complex DATA
end
t2 = t(:, index:end);

%% 2. use velocity first
% set velocity range
vel = 3.5: 0.01: 6.5;        % 理论的Sn、sSn速度是4.756；Pn是8.06419  
nvel = length(vel);

% cut window for sum
tstart = 120;      % in fact, it is tw contains phases at the reference distance !!!
tend = 140;       % CHANGE if the ref dist changes
[~, pstart] = min(abs((t1(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;

% shift and sum
datasum1 = zeros(ntw, nvel);
phasesum1 = zeros(ntw, nvel);
pwssum1 = zeros(ntw, nvel);
tuse1 = zeros(ntw, nvel);
power = 1;
figure
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
    phaseuse = phase3(pstart:pend, :);
    datasum1(:, jj) = sum(datause, 2)./nsum;      % direct sum
    phasesum1(:, jj) = abs(sum(exp(1i.*phaseuse), 2)./nsum);
    pwssum1(:, jj) = datasum1(:, jj).*(phasesum1(:, jj).^power);
    tuse1(:, jj) = t3(pstart:pend, 1);
    plot(tuse1(:, jj), datasum1(:, jj)*100 + jj, 'LineWidth',1);hold on;
end

datasum1 = datasum1./max(max(abs(datasum1)));
pwssum1 = pwssum1./max(max(abs(pwssum1)));
figure
imagesc(tuse1(:, 1), vel, datasum1');
colormap(jet);
colorbar;
caxis([-1, 1]);
xlabel('time');
ylabel('velocity');

figure
imagesc(tuse1(:, 1), vel, pwssum1');
colormap(jet);
colorbar;
caxis([-1, 1]);
xlabel('time');
ylabel('velocity');

%% 3. prepare data for shift by ray parameter P 
% NOTE: we already obtained data2 out of data, whose polarity was corrected, 
%            so we could just shift data2 according to specific phase  

% shift data2 by approximate vel. of specific phase
refvel =  5.0;                          % the reference velocity is of any specific phase 
for ii = 1:nsum
    tshift(ii) = t2(1, 1+ind4-index)-t2(1, ii)+dist(ii+index-1)./refvel;         % 8.20 km/s
    nshift(ii) =tshift(ii)/dt;      % in freq. domain, non-integer is allowed    
end
data4 = zeros(sc, nsum);
phase4 = zeros(sc, nsum);
for ii = 1:nsum
    dout = specshift(data2(:, ii), nshift(ii)-nshift(1+ind4-index));     % i(ref)=1+ind4-index, shift according to the reference trace
    data4(:, ii) = dout;      % all data has the same length as original
    pout = specshift(phase2(:, ii), nshift(ii)-nshift(1+ind4-index));
    phase4(:, ii) = pout;
end
t4= zeros(sc, nsum);
for ii = 1:nsum
    t4(:, ii) = t(:, ind4);      % all data has the same time trace, i.e. first time trace - tshift
end

figure
for ii = 1:nsum
    temp1 = data4(:, ii)./(max(abs(data4(:, ii)))) .*15.0 ;
    temp2 = (0.5*(temp1+abs(temp1)));
    plot(t4(:, ii), temp1+ double(dist(ii+index-1)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
    fill(t4(:, ii), temp2+ double(dist(ii+index-1)), 'r', 'edgealpha', 0); hold on;
%     plot(t2(figpstart:figpend, i), temp1(figpstart:figpend)+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.39 0.58 1], 'LineWidth',1); hold on;
%     fill(t2(figpstart:figpend, i), temp2(figpstart:figpend)+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
%fill([t(:, 1), fliplr(t(:, 1))], [y, fliplr(a)], 'r', 'edgealpha', 0);
end
set(gca, 'XTick', 0:5:200);
set(gca, 'xlim', [0, 200]);
set(gca, 'ylim', [390, 620]);

% % suppress the data4 from Pg/Sg to end
% suptend = 82;
% [~, aa] = min(abs(t4(:, 1)- suptend));
% x = 0:0.01:1;
% coswin = cos(pi/2.*x');
% coswin = cat(1, coswin, zeros(sc-aa+1-length(x), 1));    % design a cosine window to smoothly suppress the data to zero
% for ii=1:nsum
%     data4(aa:end, ii) = data4(aa:end, ii).*coswin;       % suppress by set the amp to 0  
% end
% % % suppress the data4 from begining until one phase
% % suptstart = 86;
% % [~, aa] = min(abs(t4(:, 1)- suptstart));
% % x = -1:0.01:0;
% % coswin = cos(pi/2.*x');
% % coswin = cat(1, zeros(aa-length(x), 1), coswin);
% % for ii=1:nsum
% %     data4(1:aa, ii) = data4(1:aa, ii).*coswin;       % suppress by set the amp to 0  
% % end
% % % ends here

figure
for ii = 1:nsum
    plot(t4(:, 1), data4(:,ii)./(max(abs(data4(:,ii)))) .*15.0 + double(dist(ii+index-1)), 'LineWidth',1);hold on;
end
set(gca,'XTick',0:2:100);

%% 4. use ray parameter P
% set parameter range
rayp = -0.06:0.001:0.06;
nrayp = length(rayp);
% time window for P parameter need to be more wide
tstart = 120; 
tend = 140;
[~, pstart] = min(abs((t4(:,1)-tstart)));
pend = pstart+round((tend-tstart)./dt);
ntw = pend-pstart+1;
datause = [];
phaseuse = [];
datasum2 = zeros(ntw, nrayp);
phasesum2 = zeros(ntw, nrayp);
pwssum2 = zeros(ntw, nrayp);
tuse2 = zeros(ntw, nrayp);
power = 1;
figure
for jj = 1:nrayp
     jj 
    for ii = 1:nsum
        tshift(ii) = dist(ii+index-1).*rayp(jj);
        nshift(ii) =tshift(ii)/dt;      % in freq. domain, non-integer is allowed    
    end
    data5 = zeros(sc, nsum);
    phase5 = zeros(sc, nsum);
    for ii = 1:nsum
        dout = specshift(data4(:, ii), nshift(ii)-nshift(1+ind4-index));     % shift according to the first trace
        data5(:, ii) = dout;      % all data has the same length as original
        pout = specshift(phase4(:, ii), nshift(ii)-nshift(1+ind4-index));
        phase5(:, ii) = pout;
    end
    t5 = t4;
    
    datause = data5(pstart:pend, :);
    datasum2(:, jj) = sum(datause, 2)./nsum;      % direct sum
    phaseuse = phase5(pstart:pend, :);
    phasesum2(:, jj) = abs(sum(exp(1i.*phaseuse), 2)./nsum);
    pwssum2(:, jj) = datasum2(:, jj).*(phasesum2(:, jj).^power);
    tuse2(:, jj) = t5(pstart:pend, 1);
    plot(tuse2(:, jj), datasum2(:, jj)*100 + jj, 'LineWidth',1);hold on;
end

% tuse2(:, 1) = tuse2(:, 1) - (545.5822-498.6431)./refvel;
datasum2 = datasum2./max(max(abs(datasum2)));
pwssum2 = pwssum2./max(max(abs(pwssum2)));
% plot datasum
figure
imagesc(tuse2(:, 1), rayp, datasum2');
colormap(jet);
colorbar;
caxis([-1, 1]);
xlabel('time');
ylabel('ray parameter');

figure
imagesc(tuse2(:, 1), rayp, pwssum2');
colormap(jet);
colorbar;
caxis([-1, 1]);
xlabel('time');
ylabel('ray parameter');
