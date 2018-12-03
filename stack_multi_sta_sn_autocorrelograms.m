%=======================================
% stacking multiple-station autocorrelograms
% reference to Zhang, Tian and Wen, GJI, 2014
% author: C. Song,   2017.3.15
%=======================================
clear ; clc ;
close all;
%% read data, real or synthetic
%cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
%datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
%datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录
%datadir = 'G:\Alxa\3test400\' ;    % 400km范围未减采数据
%datadir = 'G:\Alxa\3test500\' ;    % 500km范围未减采数据
datadir = 'G:\Alxa\3test600\' ;    % 600km范围未减采数据
sgtimedir = 'G:\Alxa\3test500_Sg_timet3\';    % 手动标记的直达Sg波到时

%fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接
fid1 = fopen(strcat(datadir,'fweight.dat')) ;
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为weight{1}  
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
% load data
[sc, sd] = size(ttest);
t = zeros(sc,sa);
data = zeros(sc,sa);
dist = zeros(sa,1);
azi = zeros(sa,1);
for i = 1:sa
    filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % 获取文件路径
    %filename = strcat(datadir,strcat(stnm(i,:),'_syn.t')) ;  % 获取文件路径
    [t(:,i), data(:,i), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(i) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
    azi(i) = SAChdr.evsta.az ;
end
% % load Sg travel time in sgtimedir
% TS = zeros(sa,1);
% for i=1:sa
%     filename = strcat(sgtimedir, strcat(stnm(i,:),'.t')) ;  % 获取文件路径
%     [~, ~, SAChdr] = fget_sac(filename) ;
%     TS(i) = SAChdr.times.t3;
% end
% TS = roundn(TS, -2);

figure (1)    % plot original data
for i = 1:sa
   % normalize
   plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
% set(gca,'XLim',[0 180]);  % 设置x轴的间隔，长度
% set(gca,'XTick',0:10:180);
% set(gca,'YLim',[30 410]);  % 设置x轴的间隔，长度
% set(gca,'YTick',30:10:410);

%% interpolate data
% dt = 0.01;
% for i=1:sa
%     t(:,i) = t(1,i): dt: t(end,i);
%     data(:,i) =  interp1(t(:,i), data(:,i), t(:,i), 'spline');
% end

%% loop for source depth to find the best one
% model parameters
n=4;     % number of layers
vs = [3.04; 3.51; 3.77; 4.756];    % S wave velocity
vp = [5.09208; 5.90612; 6.35644; 8.06419];    % P wave velocity
thick= [5.1; 20.0; 15.0; 0.0];   % thickness of each layer
nsample = 10000;     % number of samples of ray parameter P
% hmin = 5;
% hmax = 30; 
% dh = 0.1;
% nh = round((hmax-hmin)/dh+1);
% difft = zeros(sa, nh);   
tsndel = 2;    % time before direct S arrival
tsnlen = 14;     % time window length
nsnlen = round(tsnlen/dt);    % number of points of tw
thres = 5;         % threshold to control which autocorrelogram is used to stack   
sstw = 0.3;         % signal time window around zero lag
hnsstw = round(sstw/2.0/dt);     % num of points half the signal tw
nstw = 1.0;        % noise time window neighbouring
nnstw = round(nstw/dt);    % num of points the noise t
%stacksnall =zeros(2*nsnlen+1, nh);
%Edepsum = zeros(nh,1);
% start depth loop
%idep = 141;
h = 16;
% for idep = 1:nh   % depth of source
%     h = (idep-1)*dh+hmin;
    if h >= 0 && h < sum(thick(1))
        k = 1;      % layer which source locates
    elseif h >= sum(thick(1)) && h < sum(thick(1:2))
        k = 2;
    elseif h>=sum(thick(1:2)) && h<sum(thick(1:3))
        k = 3;
    else
        k = 4;
    end
    % calculate the theoratical travel time curve of each phase, like direct Sg, Sn, sSn, etc.
    % for direct S, Sg
    pS = linspace(0, 1.0/vs(k)-0.0000005, nsample)';
    if k==1 
        thicksumS = 0;
    else
        thicksumS = 0;
        for i = 1:k-1
            thicksumS = thicksumS+thick(i);
        end
    end
    XS=zeros(nsample,1);
    TS=zeros(nsample,1);
    for i = 1:nsample
        XS(i) = 0;
        hiS=zeros(k,1);
        etaS=zeros(k,1);
        for j = 1:k
            hiS(j)=thick(j);
            etaS(j) = sqrt((1.0/vs(j))^2-pS(i)^2);
            if (j == k)
                hiS(j) = h-thicksumS;
            end
            XS(i) = XS(i)+pS(i)*hiS(j)/etaS(j);
        end
        TS(i) = pS(i)*XS(i);
        for j = 1:k
            TS(i) = TS(i)+hiS(j)*etaS(j);
        end
    end

    % for Sn wave
    pSn = 1/vs(n);
    thicksumSn = 0;
    for i = 1:k
        thicksumSn = thicksumSn+thick(i);
    end
    dXSn = linspace(0, 400, nsample)';
    XSnmin = 0;
    hiSn=zeros(n-1,1);
    etaSn=zeros(n-1,1);
    for i = 1:n-1
        hiSn(i)=thick(i);
        etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
        XSnmin = XSnmin+pSn*hiSn(i)/etaSn(i);
    end
    hiSn=zeros(n-k,1);
    etaSn=zeros(n-k,1);
    for i = k:n-1
        hiSn(i)=thick(i);
        etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
        if (i == k)
            hiSn(i) = thicksumSn-h;
        end
        XSnmin = XSnmin+pSn*hiSn(i)/etaSn(i);
    end
    XSn = XSnmin + dXSn;
    TSn = pSn*XSnmin;
    hiSn=zeros(n-1,1);
    etaSn=zeros(n-1,1);
    for i = 1:n-1
        hiSn(i)=thick(i);
        etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
        TSn = TSn+hiSn(i)*etaSn(i);
    end
    hiSn=zeros(n-k,1);
    etaSn=zeros(n-k,1);
    for i = k:n-1
        hiSn(i)=thick(i);
        etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
        if (i == k)
            hiSn(i) = thicksumSn-h;
        end
        TSn = TSn+hiSn(i)*etaSn(i);
    end
    TSn = TSn + pSn*dXSn;

    % for sSn, reflected by ground surface, then like Sn 
    psSn = 1/vs(n);
    if k==1 
        thicksumsSn = 0;
    else
        thicksumsSn = 0;
        for i = 1:k-1
            thicksumsSn = thicksumsSn+thick(i);
        end
    end  
    dXsSn = linspace(0, 400, nsample)';
    XsSnmin = 0;
    hisSn=zeros(n-1,1);
    etasSn=zeros(n-1,1);
    for i = 1:n-1
        hisSn(i)=thick(i);     
        etasSn(i) = sqrt((1.0/vs(i))^2-psSn^2);    
        XsSnmin = XsSnmin+2*psSn*hisSn(i)/etasSn(i);
    end    
    hisSn=zeros(k,1);
    etasSn=zeros(k,1);
    for i = 1:k
        hisSn(i)=thick(i);
        etasSn(i) = sqrt((1.0/vs(i))^2-psSn^2);
        if (i == k)
            hisSn(i) = h-thicksumsSn;
        end        
        XsSnmin = XsSnmin+psSn*hisSn(i)/etasSn(i);
    end    
    XsSn = XsSnmin + dXsSn;
    TsSn = psSn*XsSnmin;
    hisSn=zeros(n-1,1);
    etasSn=zeros(n-1,1);
    for i = 1:n-1
        hisSn(i)=thick(i);
        etasSn(i) = sqrt((1.0/vs(i))^2-psSn^2);
        TsSn = TsSn+2*hisSn(i)*etasSn(i);
    end    
    hisSn=zeros(k,1);
    etasSn=zeros(k,1);
    for i = 1:k     
        hisSn(i)=thick(i);    
        etasSn(i) = sqrt((1.0/vs(i))^2-psSn^2);
        if (i == k)
            hisSn(i) = h-thicksumsSn;
        end
        TsSn = TsSn+hisSn(i)*etasSn(i);
    end   
    TsSn = TsSn + psSn*dXsSn;

    % interpolate  to get the  arrival time of specific distance 
    interTS = interp1(XS, TS, dist, 'spline' ); 
    interTS = roundn(interTS, -2);
    interTSn = interp1(XSn, TSn, dist, 'spline' );     % in fact, Sn, sSn appear after certain distance
    interTSn = roundn(interTSn , -2);
    interTsSn = interp1(XsSn, TsSn, dist, 'spline' );
    interTsSn = roundn(interTsSn , -2);
%     
%     % calculate the differential traveltimes of tsSn - tSn
%     if k==1 
%         thicksumsSn = 0;
%     else
%         thicksumsSn = 0;
%         for i = 1:k-1
%             thicksumsSn = thicksumsSn+thick(i);
%         end
%     end
%     hisSn=zeros(k,1);
%     etasSn=zeros(k,1);
% 
%     for i = 1:k
%         hisSn(i)=thick(i);
%         etasSn(i) = sqrt((1.0/vs(i))^2-psSn^2);
%         if (i == k)
%             hisSn(i) = h-thicksumsSn;
%         end
%         difft(:, idep) = difft(:, idep) + 2*hisSn(i)*etasSn(i);
%     end
%     difft = roundn(difft, -2);
    
    % use Sn-sSn window  to do autocorrelation
    for i=1:sa    
        if dist(i) >= 300    % when sSn is seperated from Sg
            index = i;
            break
        end        
    end   
    nuse = sa-index+1;
    
    % cut time window for Sn-sSn
    tsnb = zeros(nuse,1);    % begining of time window
    tsne = zeros(nuse,1);    % ending of time window
    indsnb =zeros(nuse,1);    % index of begining 
    tsn = zeros(nsnlen+1, nuse);     % time window
    dsn = zeros(nsnlen+1, nuse);    % data window
    figure (2)   % plot cut window
    for i =1:nuse  
        tsnb(i) = interTSn(i+index-1)-tsndel ;        % end at direct S arrival minus delay time   
        tsne(i) = tsnb(i)+tsnlen;
        [res, indsnb(i)] = min(abs(t(:, i+index-1)-tsnb(i)));
        tsn(:, i) = t(indsnb(i): indsnb(i)+nsnlen, i+index-1);   
        dsn(:, i) = data(indsnb(i): indsnb(i)+nsnlen, i+index-1);
        plot(tsn(:,i), dsn(:,i)./(max(abs(dsn(:,i)))) .*15.0 + double(dist(i+index-1)), 'linewidth', 1); hold on;
    end
    
    % use xcorr to do autocorrelation
    snlag = zeros(2*nsnlen+1, nuse);       % lag points 
    tsncorr = zeros(2*nsnlen+1, nuse);    % after correlation, length is 2*nsnlen+1
    dsncorr = zeros(2*nsnlen+1, nuse);
    figure (3)    % plot autocorrelograms
    for i = 1:nuse   
        [dsncorr(:,i), snlag(:,i)] = xcorr(dsn(:, i), 'coeff');       
        tsncorr(:,i) = snlag(:,i)*dt;
        plot(tsncorr(:,i), dsncorr(:,i)+i, 'LineWidth',1);hold on;
    end    
    set(gca,'XTick',-tsnlen: 2: tsnlen);
    
    % specify a SNR threshold, select autocorrelograms which SNR > thres
    Essn = zeros(nuse,1);
    Ensn = zeros(nuse,1);
    snrsn = zeros(nuse,1);
    sedsncorr = [];     % selected Sn-sSn data autocorrelograms to stack
    setsncorr = [];      % selected time axis
    sedistuse = [];       % selected dist
    sestnmuse = [];     % selected station
    seaziuse = [];       % selected azimuth
%     sedifftuse = [];     % selected differential time for all possible depths
    for i = 1:nuse
        % E is energy , E = sum(amp.^2)/N, zero lag is index nsnlen+1
        Essn(i) = sum(dsncorr(nsnlen+1-hnsstw : nsnlen+1+hnsstw, i).^2)./(2*hnsstw+1);    % E of signal 
        Ensn(i) = sum(dsncorr(nsnlen+1+hnsstw : nsnlen+1+hnsstw+nnstw, i).^2)/(nnstw+1);     % E of noise
        snrsn(i) = Essn(i)/Ensn(i);     % signal to nosie ratio
        if snrsn(i) > thres            
            sedsncorr = [sedsncorr dsncorr(:,i)];           
            setsncorr = [setsncorr tsncorr(:,i)];
            sedistuse = [sedistuse; dist(i+index-1)];
            seaziuse = [seaziuse; azi(i+index-1)];
            sestnmuse = [sestnmuse; stnm(i+index-1, :)];
%             sedifftuse = [sedifftuse; difft(i+index-1, idep)];    
        end        
    end    
    [~, sseu] = size(sedsncorr);
    figure (4)     % plot autocorrelograms whose SNR > threshold
    for i = 1:sseu
        plot(setsncorr(:,i), sedsncorr(:,i)*10+ sedistuse(i), 'LineWidth',1);hold on;
    end
    set(gca,'XTick', -tsnlen:2:tsnlen);
    set(gca,'YTICK', round(min(sedistuse))-5: 10: round(max(sedistuse))+5);

    % focal mechanism correction
%     strike = 266/180*pi;
%     dip = 87/180*pi;
%     rake = 5/180*pi;
%     isn = asin(vs(k)/vs(n));
%     issn = pi-isn;
%     wfact = zeros(sseu,1);
%     %figure (5)
%     for i = 1:sseu
%         % calculate radiation pattern
%         Fsn = cos(rake)*cos(dip)*cos(isn)*sin(seaziuse(i)-strike) ...
%                +cos(rake)*sin(dip)*sin(isn)*cos(2*(seaziuse(i)-strike)) ...
%                +sin(rake)*cos(2*dip)*cos(isn)*cos(seaziuse(i)-strike) ...
%                -0.5*sin(rake)*sin(2*dip)*sin(isn)*sin(2*(seaziuse(i)-strike));
%         Fssn = cos(rake)*cos(dip)*cos(issn)*sin(seaziuse(i)-strike) ...
%                 +cos(rake)*sin(dip)*sin(issn)*cos(2*(seaziuse(i)-strike)) ...
%                 +sin(rake)*cos(2*dip)*cos(issn)*cos(seaziuse(i)-strike) ...
%                 -0.5*sin(rake)*sin(2*dip)*sin(issn)*sin(2*(seaziuse(i)-strike));
%         isign = sign(Fsn)/sign(Fssn);       % sign flag 
%         if isign == 1    % same sign
%             epsilon = -1;
%         elseif isign == -1   % opposite sign
%             epsilon = 1;
%         end
%         wfact(i) = epsilon*sqrt(abs(Fsn*Fssn));
%         sedsncorr(:,i) = sedsncorr(:,i) .*wfact(i);
%         %plot(setsncorr(:,i), sedsncorr(:,i)*10+ sedistuse(i), 'LineWidth',1);hold on;
%     end
%     set(gca,'XTick', -tsnlen:2:tsnlen);
%     set(gca,'YTICK', round(min(sedistuse))-5: 10: round(max(sedistuse))+5);
    
    % stack
    % stack all
    % 1--0.5hz; 2--0.7hz; 3--0.9hz; 4--1.1hz; 5--1.3hz; 6--1.5hz
    %stacksnall = zeros(2*nsnlen+1, 6);
    stacksnall = sum(sedsncorr, 2)./sseu;
    % stack 200-300; 300-400; 400-500; 500-600
    stacksn4 = zeros(2*nsnlen+1,1);
    stacksn5 = zeros(2*nsnlen+1,1);
    stacksn6 = zeros(2*nsnlen+1,1);
    count4 = 0;    
    count5 = 0;
    count6 = 0;
    for i = 1:sseu
        if sedistuse(i) >=300 && sedistuse(i) < 400
            count4 = count4+1;
            stacksn4 = stacksn4+ sedsncorr(:,i);
        elseif sedistuse(i) >=400 && sedistuse(i) < 500
            count5 = count5+1;
            stacksn5 = stacksn5+ sedsncorr(:,i);
        else
            count6 = count6+1;
            stacksn6 = stacksn6+ sedsncorr(:,i);
        end
    end
    stacksn4 = stacksn4./count4;
    stacksn5 = stacksn5./count5;
    stacksn6 = stacksn6./count6;
    figure (6)
    if count4 >= 15, plot(setsncorr(:,1), stacksn4+2, 'LineWidth',1);hold on; end
    if count5 >= 15, plot(setsncorr(:,1), stacksn5+2.5, 'LineWidth',1);hold on; end
    if count6 >= 15, plot(setsncorr(:,1), stacksn6+3, 'LineWidth',1);hold on; end
%     plot(setsncorr(:,1), stacksn3*5+1.5, 'LineWidth',1);hold on;
%     plot(setsncorr(:,1), stacksn4*5+2, 'LineWidth',1);hold on;
%     plot(setsncorr(:,1), stacksn5*5+2.5, 'LineWidth',1);hold on;
    plot(setsncorr(:,1), stacksnall+3.5, 'LineWidth',1);hold on;
%     line([difft(1, idep), difft(1, idep)], [0, 4],'LineWidth',1);
    set(gca,'XTick',-tsnlen: 2: tsnlen);
    set(gca,'YLIM',[0 4]);
%     snall = sum(stacksnall, 2)./6;
%     figure (7)
%     for i =1:6
%         plot(setsncorr(:,1), stacksnall(:, i)*5+i/2, 'LineWidth',1);hold on;
%     end
%     plot(setsncorr(:,1), snall*5+3.5, 'LineWidth',1);hold on;
    
    % use stack energy to get best depth    
%     Edep = zeros(sseu,1);
%     for i = 1:sseu
%         ncentre = nsnlen+1+round(sedifftuse(i)./dt);    % index of theoratical sSn in each selected autocorrelogram     
%         hntw = 4;
%         Edep(i) = sum(sedsncorr(ncentre-hntw : ncentre+hntw, i).^2)./(2*hntw+1);    % energy within tw around sSn
%     end
%     Edepsum(idep) =  sum(Edep)./sseu;
% 
% %end
% [Edepmax, jbdep] = max(Edepsum);
% depbest = (jbdep-1)*dh+hmin;
% figure (7)
% for idep = 1:nh
%     plot(setsncorr(:,1), stacksnall(:, idep)*5+2.5, 'LineWidth',1);hold on;
% end
% set(gca,'XTick',-tsnlen: 2: tsnlen);
% set(gca,'YLIM',[0 8]);
% 
% hmat = hmin:dh:hmax;
% figure (8)
% plot(hmat,Edepsum);