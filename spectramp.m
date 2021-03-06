%Reads in seismograms; filters and rotates them, etc.  Can plot coherence 
%on optimal component, or amplitude spectrum on both optimal and orthogonal components,
%and their ratio, for SILB and SSIB on top of one another for specified 
%windows during one day.  
%Avrages many 1024-sample windows (at 100 sps)
format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

   %timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
   %timoffrot(1,:)=[2003 063];
   %timoffrot(1,:)=[2003 064 07 42 12 +086 +020  80 115  50];
   %timoffrot(1,:)=[2003 065 07 42 12 +086 +020  80 115  50];
%    
%    timoffrot(1,:)=[2005 254];
%    timoffrot(1,:)=[2005 255];
     %timoffrot(1,:)=[2005 228];
%    timoffrot(1,:)=[2005 256 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2005 257 07 42 12 +086 +020  80 115  50];
%   
     timoffrot(1,:)=[2004 196];
%    timoffrot(1,:)=[2004 197];
%    timoffrot(1,:)=[2004 197 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 198 07 42 12 +086 +020  80 115  50];
%    timoffrot(1,:)=[2004 199];

% % %    wins=[26045 27297; %2003.062
% % %          27916 28005;
% % %          28839 29051;
% % %          29757 29922;
% % %          38424 38789];   %wins=[26185 26275];
% % %    wins=[21800 22640; %2003.063
% % %          32300 32900;
% % %          32900 33080;
% % %          41382 41652;
% % %          41652 41832];
% % %    wins=[41473 41528]; %41503 41528]; %2003.063
% % 
  wins=[39460 40060; %2004.196
        40360 40540;
        40940 41364;
        41955 42275;
        42910 43032;
        43610 43755;
        44220 44355;
        46260 46575;
        47080 47142;
        48360 48773;
        51080 51500; 
        56230 56624;
        56845 57315;
        62500 63148;
        64935 65270;
        70420 70684;
        73940 74277;
        78885 79124;
        81520 81768;
        85830 86199];
%     wins=[1361 1387; %2004.197 (postdates 30-s)
%           5315 5620;
%           8261 8465;
%           13593 14308;
%           19920 20257;
%           25229 25625;
%           27887 27913;
%           31094 31511;
%           37815 38086;
%           39318 39597;
%           43288 43388;
%           43517 43546;
%           45911 46016;
%           46155 46259;
%           70165 70193;
%           70431 70514;
%           77370 77723;
%           77942 78043;
%           78586 78780;
%           79344 79378;
%           83567 83714;
%           84437 84554;
%           84703 84799;
%           84968 85133];
% % %   wins=[42960 42990]; %2004.196
% %     %wins=[51380 51399]; %2004.196
% % %     wins=[4245 4260; %2004.199
% % %           8968 8985;
% % %           9000 9015;
% % %           9055 9065];
% % 
% % %    wins=[14600 15320; %2005.254
% % %          21130 21400; 
% % %          27400 27580;
% % %          27580 27670;
% % %          27670 27730;
% % %          30300 30630;
% % %          30630 31080;
% % %          34500 34980;
% % %          35200 35710;
% % %          35710 35950;
% % %          49140 49920;
% % %          49920 50160;
% % %          54260 54770;
% % %          68760 69090;
% % %          73350 73800; 
% % %          79950 80370];
%      wins=[6605 7073; %2003.063
%           8830 9034;
%           10340 10456;
%           12440 12521;
%           21795 22579;
%           32350 33058;
%           41375 41692;
%           50200 50431;
%           65000 65489;
%           80425 80995];
%    wins=[14640 15270; %2005.254 
%          15740 15859; 
%          18090 18170; 
%          20035 20115; 
%          20560 20631; 
%          21155 21331; 
%          24335 24525; 
%          27430 27694;
%          30330 31065;
%          34560 34945;
%          35260 35903;
%          49140 50144;
%          54290 54719;
%          68770 69082;
%          73360 73765; 
%          79970 80340
% 	     82915 83015];
%    wins=[1000 1029; %2005.255 
% 	     5305 5800;
% 	     8060 8149;
% 	     19000 19245;
% 	     22980 23030;
%          35680 36362;
%          51030 51544;
%          54205 54359;
%          58930 59423;
%          62305 62626;
%          68200 68544;
%          75540 75638];
%      wins=[5310 5800; %2005.255
%            8064 8149;
%            19006 19245;
%            22984 23030;
%            35688 36362;
%            51033 51548;
%            54210 54359;
%            58935 59423;
%            62310 62626;
%            68205 68544;
%            75545 75638];


nd=size(timoffrot,1);

PERMSTA=['PGC'
         'LZB'];
POLSTA=['SSIB'
        'SILB'
        'KLNB'
        'TWKB'
        'MGCB'];
PERMROTS=[0 90 32 00;  %PGC , Yajun's "0" changed to 90.
          0 90 54 00]; %LZB
 POLROTS=[6 85 33 86;  %SSIB
          0 90 39 20;  %SILB
          0 90  7 00;  %KLNB
          4 70 48 00;  %MGCB
          4 75 38 00]; %TWKB
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.

yr=timoffrot(nd,1);
YEAR=int2str(yr);
jday=timoffrot(nd,2);
if jday <= 9
    JDAY=['00',int2str(jday)];
elseif jday<= 99
    JDAY=['0',int2str(jday)];
else
    JDAY=int2str(jday);
end
MO=day2month(jday,yr);
IDENTIF=[YEAR,'.',JDAY]

direc=[YEAR,'/',MO,'/'];
prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;
xero=[1.e-5 1;
      1.e5 1];
  
aaatot=zeros(7,513);
bbbtot=zeros(7,513);
ccctot=zeros(7,513);
for ista=1:size(POLSTA,1)
    hi=30;
    lo=0.01;
    sps=100;
    STAEdat=[prename,'.',POLSTA(ista,:),'..HHE.D.SAC']; %BHE for permstas.
    STANdat=[prename,'.',POLSTA(ista,:),'..HHN.D.SAC'];
    [STAE,HdrDataSTA,tnuSTA,pobjSTA,timsSTA]=readsac(STAEdat,0,'l');
    [STAN,~,~,~,~]=readsac(STANdat,0,'l');
    tracelen=length(STAE);
    %cosine taper before filtering:
    x=(0:pi/200:pi/2-pi/200)';
    STAE(1:100)=sin(x).*STAE(1:100);
    STAN(1:100)=sin(x).*STAN(1:100);
    x=flipud(x);
    STAE(tracelen-99:tracelen)=sin(x).*STAE(tracelen-99:tracelen);
    STAN(tracelen-99:tracelen)=sin(x).*STAN(tracelen-99:tracelen);
    %Filter data:
    npo=2;
    npa=1;
    if yr==2003 && jday<213
        [STAEf]=20.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
        [STANf]=20.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); %times 5 until Aug 2003 5.0*
    else
        [STAEf]=4.0e-3*bandpass(STAE,sps,lo,hi,npo,npa,'butter'); 
        [STANf]=4.0e-3*bandpass(STAN,sps,lo,hi,npo,npa,'butter'); 
    end
    STA=STAEf+1i*STANf;
    STAfastslow=STA*exp(-1i*POLROTS(ista,2));
    STAslow=real(STAfastslow);
    STAfast=imag(STAfastslow);
    len=length(STA);
    STAslow(20:len-20)=STAslow(20+POLROTS(ista,1):len-20+POLROTS(ista,1));
    
    STAsplitcorrected=(STAslow+1i*STAfast)*exp(1i*POLROTS(ista,2));
    STAsplitcorrectedrot=STAsplitcorrected*exp(-1i*POLROTS(ista,3));
    STAscrot=STAsplitcorrectedrot;

%     SSIBEfd = resample(SSIBEf,2,5);
%     SSIBNfd = resample(SSIBNf,2,5);
%     SILBEfd = resample(SILBEf,2,5);
%     SILBNfd = resample(SILBNf,2,5);

    STAsoff=POLROTS(ista,4)
    STAtoff=round(STAsoff/40); %really 40.  That's how it's listed.

    if STAtoff > -1
        STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
        STAscrot(tracelen-STAsoff+1:tracelen)=0;
    else
        STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
        STAscrot(1:-STAsoff)=0;
    end
    STAopt=real(STAscrot);
    STAort=imag(STAscrot);
    STAoptfilt=bandpass(STAopt,sps,1.5,6,npo,npa,'butter');
    STAoptenv=abs(hilbert(STAoptfilt));
    
    swin=1024; %for small window
    itot=0;
    for n=1:size(wins,1)
        winstart=wins(n,1)*sps;
        winend=wins(n,2)*sps;
        nswins=2*floor((winend-winstart)/swin);
        for iswin=1:nswins
            itot=itot+1;
            istart=winstart+(iswin-1)*(swin/2);
            iend=istart-1+swin;
            [STAoptxx(itot,1:swin/2+1) f1] = periodogram(STAopt(istart:iend),[],[],sps);
            [STAortxx(itot,1:swin/2+1) f1] = periodogram(STAort(istart:iend),[],[],sps);
            STAoptxx(itot,swin/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
            STAortxx(itot,swin/2+2) = mean(STAoptenv(istart:iend));  %Last row gets estimate of envelope amplitude (optimal orientation)
            %plot(log10(f1),log10(STAoptxx))
            %loglog(f1,STAoptxx)
            %drawnow
            %[STAortxx f1] = periodogram(STAort(istart:iend),[],[],sps);
            %loglog(f1,STAortxx,'r')
            %loglog(f1,STAoptxx./STAortxx,'g')
            %loglog(xero(:,1),xero(:,2),'k')
            %loglog(xero(:,1),10*xero(:,2),'k')
            %loglog(xero(:,1),1.e-4*xero(:,2),'k')
            
        end
%         title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ',int2str(wins(n,1)),'-',int2str(wins(n,2)),' sec'])
%         set(h,'PaperPosition',[0.25 5 8 5])
%         print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'_',int2str(wins(n,1)),'-',int2str(wins(n,2)),'.eps'])
    end
    a=sortrows(STAoptxx,514); %Each row is a periodogram; first row smallest amplitude; last row largest.
    aa=cumsum(a(:,1:513));  %cumsum cumsums the columns
    b=sortrows(STAortxx,514);
    bb=cumsum(b(:,1:513));
    h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
    subplot(3,1,1,'align')
%     aaa(1,:)=aa(64,:)/64;  %The first 64 rows (aa is a cumsum)
%     for i=2:floor(itot/64)
%         aaa(i,:)=(aa(i*64,:)-aa((i-1)*64,:))/64;
%     end
%     aaa(i+1,:)=(aa(end,:)-aa(i*64,:))/(itot-i*64);
    ends=110;
    divs=[ends/2 round(0.2*(itot-ends)) round(0.4*(itot-ends)) round(0.6*(itot-ends)) round(0.8*(itot-ends)) itot-ends/2 itot];
    aaa(1,:)=aa(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
    for i=2:length(divs)
        aaa(i,:)=(aa(divs(i),:)-aa(divs(i-1),:))/(divs(i)-divs(i-1));
    end
    if ~strcmp(POLSTA(ista,:),'MGCB')
        aaatot=aaa+aaatot;
    end
    loglog(f1,aaa)
    hold on
    load SILBaaa.mat
    loglog(f1,aaanoise,'k--')
    title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  opt'])
    xlim([0.1 30])
    ylim([1.e-7 1])
    subplot(3,1,2,'align')
%     bbb(1,:)=bb(64,:)/64;
%     for i=2:floor(itot/64)
%         bbb(i,:)=(bb(i*64,:)-bb((i-1)*64,:))/64;
%     end
%     bbb(i+1,:)=(bb(end,:)-bb(i*64,:))/(itot-i*64);
    bbb(1,:)=bb(divs(1),:)/divs(1);  %The first 64 rows (aa is a cumsum)
    for i=2:length(divs)
        bbb(i,:)=(bb(divs(i),:)-bb(divs(i-1),:))/(divs(i)-divs(i-1));
    end
    if ~strcmp(POLSTA(ista,:),'MGCB')
        bbbtot=bbb+bbbtot;
    end
    loglog(f1,bbb)
    hold on
    load SILBbbb.mat
    loglog(f1,bbbnoise,'k--')
    title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  ort'])
    xlim([0.1 30])
    ylim([1.e-7 1])
    subplot(3,1,3,'align')
    ccc=aaa./bbb;
    if ~strcmp(POLSTA(ista,:),'MGCB')
        ccctot=ccc+ccctot;
    end
    loglog(f1,ccc)
    hold on
    title([POLSTA(ista,:),'  ',YEAR,'.',JDAY,'  opt/ort'])
    xlim([0.1 30])
    ylim([1.e-1 11])
    loglog([1.e-5 1.e5],[1 1],'k--')
    set(h,'PaperPosition',[0.25 0.25 8 10])
    print(h,'-depsc',['spectr',POLSTA(ista,:),YEAR,'.',JDAY,'.eps'])
    clear STAoptxx
    clear STAortxx
end

h=figure('Position',[wid/10 hite/9 1.*wid hite/1.2]);
% allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' ',POLSTA(4,:),' ',POLSTA(5,:),' '];
% nsta=5;
% allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' ',POLSTA(5,:),' '];
% nsta=4;
allsta=[POLSTA(1,:),' ',POLSTA(2,:),' ',POLSTA(3,:),' ',POLSTA(4,:),' '];
nsta=4;
subplot(3,1,1,'align')
loglog(f1,aaatot/nsta)
hold on
title([allsta,'  ',YEAR,'.',JDAY,'  opt'])
xlim([0.1 30])
ylim([1.e-7 1])
subplot(3,1,2,'align')
loglog(f1,bbbtot/nsta)
title([allsta,'  ',YEAR,'.',JDAY,'  ort'])
xlim([0.1 30])
ylim([1.e-7 1])
subplot(3,1,3,'align')
loglog(f1,ccctot/nsta)  %size(POLSTA,1)
hold on
loglog([1.e-5 1.e5],[1 1],'k--')
title([allsta,'  ',YEAR,'.',JDAY,'  opt/ort'])
xlim([0.1 30])
ylim([1.e-1 11])
set(h,'PaperPosition',[0.25 0.25 8 10])
print(h,'-depsc',['spectr',YEAR,'.',JDAY,'.eps'])

