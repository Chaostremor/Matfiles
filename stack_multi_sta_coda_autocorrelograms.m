%=======================================
% stacking multiple-station autocorrelograms use S coda wave 
% reference to Zhang, Tian and Wen, GJI, 2014
% author: C. Song,   2017.3.15
%=======================================
clear ; clc ; close all;
%% read data, real or synthetic
%cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
%datadir = 'G:\Alxa\400_15km_real\' ;   % 数据所在目录
%datadir = 'G:\Alxa\400_15km_syn\' ;    % 理论所在目录
%datadir = 'G:\Alxa\3test400\' ;    % 400km范围未减采数据
datadir = 'G:\Alxa\3test500\' ;    % 500km范围未减采数据
sgtimedir = 'G:\Alxa\3test500_Sg_timet3\';    % 手动标记的直达Sg波到时
%fid1 = fopen(strcat(datadir,'fweight_select.dat')) ;      % strcat用于字符串连接
fid1 = fopen(strcat(datadir,'fweight.dat')) ;
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为weight{1}  
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

%% interpolate data
% dt = 0.01;
% for i=1:sa
%     t(:,i) = t(1,i): dt: t(end,i);
%     data(:,i) =  interp1(t(:,i), data(:,i), t(:,i), 'spline');
% end

%%  calculate the theoratical travel time curve of each phase, like direct Sg, Sn, sSn, etc.
% model parameters
n=4;     % number of layers
vs = [3.04; 3.51; 3.77; 4.756];    % S wave velocity
vp = [5.09208; 5.90612; 6.35644; 8.06419];    % P wave velocity
thick= [5.1; 20.0; 15.0; 0.0];   % thickness of each layer
h = 16.0;     % depth of source
if h >= 0 && h < sum(thick(1))
    k = 1;      % layer which source locates
elseif h >= sum(thick(1)) && h < sum(thick(1:2))
    k = 2;
elseif h>=sum(thick(1:2)) && h<sum(thick(1:3))
    k = 3;
else
    k = 4;       % in fact, if k==4, the whole algorithm is not useful,
end              % we talk about sources above Moho
nsample = 10000;     % number of samples of ray parameter P
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

%% interpolate  to get the  arrival time of specific distance 
interTS = interp1(XS, TS, dist, 'spline' );
interTS = roundn(interTS, -2);
interTSn = interp1(XSn, TSn, dist, 'spline' );
interTSn = roundn(interTSn , -2);
interTsSn = interp1(XsSn, TsSn, dist, 'spline' );
interTsSn = roundn(interTsSn , -2);

%% use coda wave data to do autocorrelation
% coda wave is selected as 30-60s after the first theoratical arrival time of S wave,
% which means Sg when dist is small, Sn when dist is large, about 200 km 
%for j = 1:6
%tcodadel = 40+(j-1)*5
tcodadel = 60;    % delay time from  S arrival
tcodalen = 30;     % time window length
ncodalen = round(tcodalen/dt);    % number of points of tw
tcodab = zeros(sa,1);    % begining of time window
tcodae = zeros(sa,1);    % ending of time window
indcodab =zeros(sa,1);    % index of begining 
tcoda = zeros(ncodalen+1,sa);     % time window
dcoda = zeros(ncodalen+1,sa);    % data window
for i =1:sa
    %tcodab(i) = min([interTS(i) interTSn(i)])+tcodadel ;     % start from first S arrival plus delay time
    tcodab(i) = interTS(i)+tcodadel ;        % start from direct S arrival plus delay time
    tcodae(i) = tcodab(i)+tcodalen;
    [res, indcodab(i)] = min(abs(t(:,i)-tcodab(i)));
    tcoda(:, i) = t(indcodab(i):indcodab(i)+ncodalen, i);
    dcoda(:, i) = data(indcodab(i):indcodab(i)+ncodalen, i);
end
% use xcorr to do autocorrelation
codalag = zeros(2*ncodalen+1,sa);       % lag points 
tcodacorr = zeros(2*ncodalen+1,sa);    % after correlation, length be 2*ncodalen+1
dcodacorr = zeros(2*ncodalen+1,sa);
figure (2)
for i = 1:sa
    [dcodacorr(:,i), codalag(:,i)] = xcorr(dcoda(:,i),'coeff');
    tcodacorr(:,i) = codalag(:,i)*dt;
    plot(tcodacorr(:,i), dcodacorr(:,i)+i, 'LineWidth',1);hold on;
end
set(gca,'XTick',-tcodalen: 2: tcodalen);
% specify a SNR threshold, select autocorrelograms which SNR > thres
thres = 5;         % threshold
sctw = 0.3;         % signal time window of coda
hnsctw = round(sctw/2.0/dt);     % num of points half the signal tw of coda
nctw = 1.0;        % noise time window of coda
nnctw = round(nctw/dt);    % num of points the noise tw of coda
Esco = zeros(sa,1);
Enco = zeros(sa,1);
snrco = zeros(sa,1);
sedcocorr = [];     % selected coda data autocorrelograms to stack
setcocorr = [];      % selected coda time axis
sedist = [];       % selected dist
sestnm = [];     % selected station 
for i = 1:sa
    % E is energy , E = sum(amp.^2)/N
    Esco(i) = sum(dcodacorr(ncodalen+1-hnsctw : ncodalen+1+hnsctw, i).^2)./(2*hnsctw+1);    % E is energy, zero lag is index ncodalen+1
    Enco(i) = sum(dcodacorr(ncodalen+2+hnsctw : ncodalen+2+hnsctw+nnctw, i).^2)/(nnctw+1);
    snrco(i) = Esco(i)/Enco(i);
    if snrco(i) > thres
        sedcocorr = [sedcocorr dcodacorr(:,i)];
        setcocorr = [setcocorr tcodacorr(:,i)];
        sedist = [sedist; dist(i)];
        sestnm = [sestnm; stnm(i,:)];
    end
end
[~, sse] = size(sedcocorr);
figure (3)
for i = 1:sse
    plot(setcocorr(:,i), sedcocorr(:,i)*10+ sedist(i), 'LineWidth',1);hold on;
end
set(gca,'XTick',-tcodalen: 2: tcodalen);
set(gca,'YTICK', round(min(sedist))-5: 10: round(max(sedist))+5);
% stack the autocorrelograms 
% stack all;
stackcoall = sum(sedcocorr, 2)./sse;
% stack 0-100km; 100-200; 200-300; 300-400; 400-500;
stackco1 = zeros(2*ncodalen+1,1);
stackco2 = zeros(2*ncodalen+1,1);
stackco3 = zeros(2*ncodalen+1,1);
stackco4 = zeros(2*ncodalen+1,1);
stackco5 = zeros(2*ncodalen+1,1);
count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
count5 = 0;
for i = 1:sse
    if sedist(i) >0 && sedist(i) <= 100
        count1 = count1+1;
        stackco1 = stackco1+ sedcocorr(:,i);
    elseif sedist(i) >100 && sedist(i) <= 200
        count2 = count2+1;
        stackco2 = stackco2+ sedcocorr(:,i);
    elseif sedist(i) >200 && sedist(i) <= 300
        count3 = count3+1;
        stackco3 = stackco3+ sedcocorr(:,i);
    elseif sedist(i) >300 && sedist(i) <= 400
        count4 = count4+1;
        stackco4 = stackco4+ sedcocorr(:,i);
    else
        count5 = count5+1;
        stackco5 = stackco5+ sedcocorr(:,i);
    end
end
stackco1 = stackco1./count1;
stackco2 = stackco2./count2;
stackco3 = stackco3./count3;
stackco4 = stackco4./count4;
stackco5 = stackco5./count5;
% if count1 < 15
%     stackco1 = 0;
% end
% if count2 < 15
%     stackco2 = 0;
% end
% if count3 < 15
%     stackco3 = 0;
% end
% if count4 < 15
%     stackco4 = 0;
% end
% if count5 < 15
%     stackco5 = 0;
% end

[max1, imax1] = max(stackcoall(ncodalen+401:end));
[min1, imin1] = min(stackcoall(ncodalen+401:end));
tmax = setcocorr(ncodalen+400+imax1, 1);
tmin = setcocorr(ncodalen+400+imin1, 1);

nn = 30;
Esmax = sum(stackcoall(ncodalen+400+imax1-nn : ncodalen+400+imax1+nn).^2)./(2*nn+1);
Enmax = (sum(stackcoall(ncodalen+400+imax1-3*nn-1: ncodalen+400+imax1-nn-1).^2)./(2*nn+1)+sum(stackcoall(ncodalen+400+imax1+nn+1: ncodalen+400+imax1+3*nn+1).^2)./(2*nn+1))/2;
SNRmax = Esmax/Enmax;
Esmin = sum(stackcoall(ncodalen+400+imin1-nn : ncodalen+400+imin1+nn).^2)./(2*nn+1);
Enmin = (sum(stackcoall(ncodalen+400+imin1-3*nn-1: ncodalen+400+imin1-nn-1).^2)./(2*nn+1)+sum(stackcoall(ncodalen+400+imin1+nn+1: ncodalen+400+imin1+3*nn+1).^2)./(2*nn+1))/2;
SNRmin = Esmin/Enmin;

figure (4)
if count1 >= 15, plot(setcocorr(:,1), stackco1+0.5, 'LineWidth',1);hold on; end
if count2 >= 15, plot(setcocorr(:,1), stackco2+1, 'LineWidth',1);hold on; end
if count3 >= 15, plot(setcocorr(:,1), stackco3+1.5, 'LineWidth',1);hold on; end
if count4 >= 15, plot(setcocorr(:,1), stackco4+2, 'LineWidth',1);hold on; end
if count5 >= 15, plot(setcocorr(:,1), stackco5+2.5, 'LineWidth',1);hold on; end
plot(setcocorr(:,1), stackcoall+3, 'LineWidth',1);hold on;
line([tmax, tmax], [0, 4],'LineWidth',1);hold on;
line([tmin, tmin], [0, 4],'LineWidth',1);
set(gca,'XTick',-tcodalen: 2: tcodalen);
set(gca,'YLIM',[0 4]);
%end

% figure (5)
% for i = 1:6
%     plot(setcocorr(:,1), stackcoall(:, i)+i/5, 'LineWidth',1);hold on;
% end
% plot(setcocorr(:,1), sum(stackcoall, 2)/6+1.4, 'LineWidth',1);
% set(gca,'XTick',-tcodalen: 2: tcodalen);
% set(gca,'YLIM',[0 1.6]);