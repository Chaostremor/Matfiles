format short e
clear all
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');
scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

%timoffrot(1,:)=[2003 062 07 42 12 +086 +020  80 120  50]; %120 for 2-8 Hz
%timoffrot(1,:)=[2003 063 07 42 12 +086 +020  80 115  50];
%timoffrot(3,:)=[2003 064 07 42 12 +086 +020  80 115  50];
timoffrot= [2003 063];
timoffrot= [2005 255];
% timoffrot=[2003 062;
%            2003 063;
%            2003 064;
%            2003 068;
%            2004 196;
%            2004 197;
%            2004 198;
%            2004 199;
%            2004 200;
%            2005 254;
%            2005 255;
%            2005 256;
%            2005 260];

Bost='030304';
Bost='050912';
% Bost=['030303'
%       '030304'
%       '030305'
%       '030309'
%       '040714'
%       '040715'
%       '040716'
%       '040717'
%       '040718'
%       '050911'
%       '050912'
%       '050913'
%       '050917'];

PERMSTA=['PGC  '
         'LZB  '
         'LZBz '];
POLSTA =['SSIB '
         'SILB '
         'KLNB '
         'MGCB '
         'TWKB '
         'SILBz'
         'SSIBz'];
fam='002';
% 002
PERMROTS=[0 90 32 0;  %PGC , Yajun's "0" changed to 90.
          0 90 54 9;  %LZB
          0  0  0 9]; %LZBz
POLROTS=[6 85 33 86;  %SSIB from Yajun
         0 90 39 20;  %SILB
         0 90  7 -4;  %KLNB
         4 70 48 -26; %MGCB
         4 75 38 -5;  %TWKB
         0  0  0 20;  %SILBz
         0  0  0 86]; %SSIBz
% POLROTS=[6 85 33 86;  %SSIB from Yajun
%          0 90 39 20;  %SILB
%          0 90  7 -5;  %KLNB
%          4 70 48 -4;  %MGCB
%          4 75 38 -27]; %TWKB
% Column order is:
%   1-Bostock family ID# 
%   2-fast/slow time offset (at 40 sps)
%   3-slow direction 
%   4-fast/slow cc coefficient 
%   5-corrected polarization angle
%   6-uncorrected polarization angle.  Angles are c-clockwise from East
% Station order is 1.PGC 2.LZB 3.VGZ 4.SSIB 5.SILB 6.TWKB 7.MGCB 8.TWBB 9.TSJB 10.KLNB
% 2.0000  0.0000  0.0000  0.0000  32.0000 22.8366
% 2.0000  0.0000  0.0000  0.0000  55.0000 17.7860
% 2.0000  0.0000  0.0000  0.0000  96.0000 2.5716
% 2.0000  -6.0000 85.0000 -0.9470 33.0000 41.9115
% 2.0000  0.0000  0.0000  0.0000  40.0000 23.3330
% 2.0000  -4.0000 75.0000 -0.9393 38.0000 36.5957
% 2.0000  -4.0000 70.0000 -0.9242 46.0000 114.1442
% 2.0000  0.0000  0.0000  0.0000  0.0000  0.0000
% 2.0000  0.0000  0.0000  0.0000  0.0000  0.0000
% 2.0000  0.0000  0.0000  0.0000  8.0000  9.0570

% 047
%  POLROTS=[4 80 29 86;  %SSIB
%           0 90 41 20;  %SILB
%           0 90 00 00;  %KLNB
%           2 85 49 00;  %MGCB
%           4 65 24 00]; %TWKB
% 068
%  POLROTS=[2 50 25 00; %KLNB
%           4 60 14 00; %TWKB
%           5 65 41 00]; %MGCB
stas=['PGC  '
      'SSIB '
      'SILB '];
% stas=['LZB  '
%       'TWKB '
%       'MGCB '];
nsta=size(stas,1);
PERMROTS(:,2:3)=pi*PERMROTS(:,2:3)/180.;
POLROTS(:,2:3)=pi*POLROTS(:,2:3)/180.;
%POLROTS(:,1)=round(POLROTS(:,1)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
%POLROTS(:,4)=round(POLROTS(:,4)*(100/40)); %Polaris stations at 100 sps; table has split time at 40 sps.
sps=40;
PERMROTS(:,1)=round(PERMROTS(:,1)*(sps/40));  %For some reason PERMROTS(:,4) is still adjusted in readperms.m
POLROTS(:,1)=round(POLROTS(:,1)*(sps/40));
tempwinlen=60*sps;
stack=zeros(nsta,tempwinlen);
stackort=zeros(nsta,tempwinlen);

%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensecCC=2;
winlensec=8;
winoffsec=1; %I think readpols etc. need this for counting glitches still.
winlenCC=winlensecCC*sps;
winlen=winlensec*sps;
winoff=winoffsec*sps;
tracelen=86400*sps; %one day of data at 40 sps
winbig=2*(tracelen/2-(2*sps)); %ignore 2 seconds at each end of day
timbig=winbig/(2*sps); %half that time, in seconds
igstart=floor(tracelen/2-winbig/2)+1; %start counting seis data from here
nwin=floor((winbig-winlen)/winoff);
hi=6.5;
lo=1.25;
npo=2;
npa=2;
mshift=round(10*sps/40); %maximum shift for the x-correlations.
loopoffmax=1.5*sps/40; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.
xcmaxAVEnmin=0.45; %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?

%cycle over each day:
for nd=1:length(timoffrot(:,1))
    close all
    clear AmpComp
%     %Bostock's detections:
%     bostocks=load(bostname(nd,:));
%     bostsec=3600*(bostocks(:,3)-1)+bostocks(:,4)-22.675; %22.675 b.c. 002_ catalogs already "corrected" for PGC start; from tempseps.f
%     bostsamp=round(bostsec*40);
    %Which days of data to read?
    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    MO=day2month(jday,year);
    IDENTIF=[YEAR,'.',JDAY,'.',fam];
    direc=[YEAR,'/',MO,'/'];
    prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];
    
    %Read times of Bostock's detections:
    % 2  30303  6 1158.350 1.427 13
    % 2  30303  7  650.050 1.468 14
    % 2  30303  7  781.900 1.339 16
    % 2  30303  7  812.300 1.516 13
    BostCat=load(['/data2/arubin/CNDC/BOSTOCK/NEW/',Bost(nd,:)]);
    timswin=3600.*(BostCat(:,3)-1)+BostCat(:,4)+0.076; %+0.076 to correct from 22.675s to 22.751s arrival (2-pass 1.5-7 Hz) 
    nwin=length(timswin);
    iarrive=round(sps*timswin);
    istartCC=iarrive-round(winlenCC/2); %Run the CC over a shorter time window than plotted.
    iendCC=istartCC+winlenCC-1;
    istart=iarrive-round(3*winlen/4); %can off-center this to pick up P rather than S coda.
    iend=istart+winlen-1;
    iccwinoff=istartCC(1)-istart(1);

    %Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
    %Get timsSTA from the permanent stations (last one over-writes):
    STAopt=zeros(nsta,tracelen);
    STAort=STAopt;
    STAnzeros=zeros(nsta,nwin);
    for ista=1:nsta
        found=0;
        [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [opt,ort,nzeros,timsSTAperm]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort,nzeros,timsSTApol]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        found=found
        STAopt(ista,:)=opt;
        %STAort(ista,:)=ort;tim
        STAnzeros(ista,:)=nzeros;
    end
    if sps==40
        timsSTA=timsSTAperm;
    else
        timsSTA=timsSTApol;
    end
    clear ort
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Autocorrelation of stations.  Those that end in "sq" are the running
    %   cumulative sum, to be used later by differncing the window edpoints.
    %   (Used to be PGCauto, PGC2, SSIBauto, SSIB2, etc.)
    %   Station to itself is in a 3 x tracelen array
    %   Cross-station measurements are in their own linear array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    STAauto=STAopt.*STAopt;
    STAsq=cumsum(STAauto,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Cross-correlation between stations, with small offsets up to +/- mshift.
    %  First index is pointwise multiplication of traces; second is shifting offset.
    %  lenx is shorter than tracelen by mshift at each end (see notebook sketch)
    %  For STA12 and PGSI, SSI and SIL are shifted relative to PGC, by 1 each time through loop.
    %  For SISS, SSI is shifted relative to SILB.
    %  cumsumSTA12 etc. are the running cumulative sum of the x-correlation.
    %  PGSSx becomes STA12x, PGSI -> STA13, SISS -> STA32
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lenx=tracelen-2*mshift;
    STA12x=zeros(lenx, 2*mshift+1);
    STA13x=zeros(lenx, 2*mshift+1);
    STA32x=zeros(lenx, 2*mshift+1);
    for n=-mshift:mshift;
        STA12x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
        STA13x(:,n+mshift+1)=STAopt(1,1+mshift:tracelen-mshift).* ...
            STAopt(3,1+mshift-n:tracelen-mshift-n);
        STA32x(:,n+mshift+1)=STAopt(3,1+mshift:tracelen-mshift).* ...
            STAopt(2,1+mshift-n:tracelen-mshift-n);
    end
    cumsumSTA12=cumsum(STA12x);  %prev cumsumPGSS
    cumsumSTA13=cumsum(STA13x);  %prev cumsumPGSI
    cumsumSTA32=cumsum(STA32x);  %prev cumsumSISS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  "winbig" is now the whole day, minus 2 sec at each end (apparently).
    %  "timbig" is the time of half that.
    %  igstart is the index of the starting sample.
    %  winlen and winoff refer to the small windows.
    %  timswin refers to the central times of those small windows.
    %  sumsPGSS (etc.) is the cross-correlation sum over the window.  The first
    %    index refers to the window number and the second the shift over +/-mshift.
    %  Normalized x-correlation:
    %    For PGSS and PGSI, for a given window PGC does not shift but SSI and 
    %    SIL do.  So can compute sumsPGC2 (from the running cum. sum PGC2) just
    %    once for each window.  Same for sumsSILB2b.  But for the stations that
    %    shift, SSI and SIL (for PGC) and SSI (for SIL), must compute sumsSSIB2 
    %    and sumsSILB2 upon each shift (actually, this is is easy book-keeping
    %    but not efficient).  Again, the first index refers to the window
    %    number and the second the shift over +/-mshift.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sumsSTA12=zeros(nwin,2*mshift+1);
    sumsSTA13=zeros(nwin,2*mshift+1);
    sumsSTA32=zeros(nwin,2*mshift+1);
    sumsSTA1sq=zeros(nwin,2*mshift+1);
    sumsSTA2sq=zeros(nwin,2*mshift+1);
    sumsSTA3sq=zeros(nwin,2*mshift+1);
    sumsSTA3Bsq=zeros(nwin,2*mshift+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  sumsPGSS is shorter than sumsPGC2 by 2*mshift.  This is why sumsPGC2 etc
    %  is shifted by +mshift.  cumsumPGSS(1,:)=cumsum(PGSSx)(1,:) starts mshift
    %  to the right of the first data sample.  igstart is how many to the right
    %  of that.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n=1:nwin;
%         istart=igstart+(n-1)*winoff;
%         iend=istart+winlen;
%         timswin(n)=timsSTA(istart+winlen/2); 
        sumsSTA12(n,:)=cumsumSTA12(iendCC(n),:)-cumsumSTA12(istartCC(n)-1,:); 
        sumsSTA13(n,:)=cumsumSTA13(iendCC(n),:)-cumsumSTA13(istartCC(n)-1,:);
        sumsSTA32(n,:)=cumsumSTA32(iendCC(n),:)-cumsumSTA32(istartCC(n)-1,:);
        sumsSTA1sq(n,:)=STAsq(1,iendCC(n)+mshift)-STAsq(1,istartCC(n)+mshift-1);  %PGC2 is cumsummed. Yes, +mshift.
        sumsSTA3Bsq(n,:)=STAsq(3,iendCC(n)+mshift)-STAsq(3,istartCC(n)+mshift-1); %Similar, for the SILB-SSIB connection.
        for m=-mshift:mshift;
            sumsSTA2sq(n,m+mshift+1)=STAsq(2,iendCC(n)+mshift-m)-STAsq(2,istartCC(n)+mshift-1-m); %+m??? (yes).
            sumsSTA3sq(n,m+mshift+1)=STAsq(3,iendCC(n)+mshift-m)-STAsq(3,istartCC(n)+mshift-1-m);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Denominator for the normalization.  A 2D array, nwin by 2*mshift+1.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %An attempt to bypass glitches in data.  Min value of good data typically ~10^{-2}
    glitches=1.e-7;
    sumsSTA1sq=max(sumsSTA1sq,glitches);
    sumsSTA2sq=max(sumsSTA2sq,glitches);
    sumsSTA3sq=max(sumsSTA3sq,glitches);
    %
    denomSTA12n=realsqrt(sumsSTA1sq.*sumsSTA2sq);
    denomSTA13n=realsqrt(sumsSTA1sq.*sumsSTA3sq);
    denomSTA32n=realsqrt(sumsSTA3Bsq.*sumsSTA2sq);
    %
    sumsSTA12n=sumsSTA12./denomSTA12n;
    sumsSTA13n=sumsSTA13./denomSTA13n;
    sumsSTA32n=sumsSTA32./denomSTA32n;
    [xcmaxSTA12n,imaxSTA12]=max(sumsSTA12n,[],2); %Integer-offset max cross-correlation
    [xcmaxSTA13n,imaxSTA13]=max(sumsSTA13n,[],2);
    [xcmaxSTA32n,imaxSTA32]=max(sumsSTA32n,[],2);
    %Parabolic fit:
    [xmaxSTA12n,ymaxSTA12n,aSTA12]=parabol(nwin,mshift,sumsSTA12n,imaxSTA12); %Interpolated max cross-correlation
    [xmaxSTA13n,ymaxSTA13n,aSTA13]=parabol(nwin,mshift,sumsSTA13n,imaxSTA13);
    [xmaxSTA32n,ymaxSTA32n,aSTA32]=parabol(nwin,mshift,sumsSTA32n,imaxSTA32);
    
    ix=sub2ind(size(denomSTA12n),(1:nwin)',imaxSTA12); %Find the linear index of the largest denominator
    ampSTA12=sqrt(denomSTA12n(ix)); %This makes amplitude linear rather than quadratic with counts.
    ampSTA1sq=sumsSTA1sq(ix); %by construction PGC2 is the same for all shifts  % sumsPGC2 becomes sumsSTA1sq
    ampSTA2sq=sumsSTA2sq(ix); % sumsSSIB2 becomes sumsSTA2sq
    ix=sub2ind(size(denomSTA13n),(1:nwin)',imaxSTA13);
    ampSTA13=sqrt(denomSTA13n(ix));
    ampSTA3sq=sumsSTA3sq(ix);
    ix=sub2ind(size(denomSTA32n),(1:nwin)',imaxSTA32);
    ampSTA32=sqrt(denomSTA32n(ix));
    AmpComp(1:4)=0;
    AmpComp(5:nwin)=((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))- ...
                    (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4)))./ ...
                    ((ampSTA1sq(5:nwin)+ampSTA2sq(5:nwin)+ampSTA3sq(5:nwin))+ ...
                    (ampSTA1sq(1:nwin-4)+ampSTA2sq(1:nwin-4)+ampSTA3sq(1:nwin-4))) ;
    %Center them
    imaxSTA12cent=imaxSTA12-mshift-1;  % "cent" is "centered"; imaxSTA12 is original 1:2*mshift+1
    imaxSTA13cent=imaxSTA13-mshift-1;
    imaxSTA32cent=imaxSTA32-mshift-1;
    iloopoff=imaxSTA13cent-imaxSTA12cent+imaxSTA32cent; %How well does the integer loop close?
    xmaxSTA12n=xmaxSTA12n-mshift-1;
    xmaxSTA13n=xmaxSTA13n-mshift-1;
    xmaxSTA32n=xmaxSTA32n-mshift-1;
    loopoff=xmaxSTA13n-xmaxSTA12n+xmaxSTA32n; %How well does the interpolated loop close?
    xcmaxAVEn=(xcmaxSTA12n+xcmaxSTA13n+xcmaxSTA32n)/3;
    % xcnshifts=cputime-t
    % t=cputime;
    ampmax=max([ampSTA12; ampSTA13; ampSTA32]);
    medxcmaxAVEn=median(xcmaxAVEn)
    xmaxSTA12ntmp=xmaxSTA12n;
    xmaxSTA13ntmp=xmaxSTA13n;
    xmaxSTA32ntmp=xmaxSTA32n;

    % iup=4;
    nin=0;
    zerosallowed=20*winlen/(4*sps);
    for n=1:nwin
%         if xcmaxAVEn(n)<xcmaxAVEnmin || abs(xmaxSTA13n(n)-xmaxSTA12n(n)+xmaxSTA32n(n))>loopoffmax ...
%                 || isequal(abs(imaxSTA12cent(n)),mshift) || isequal(abs(imaxSTA13cent(n)),mshift) ...
%                 || isequal(abs(imaxSTA32cent(n)),mshift) || max(STAnzeros(:,n))>zerosallowed        
%             xmaxSTA12ntmp(n)=mshift+1; xmaxSTA13ntmp(n)=mshift+1; xmaxSTA32ntmp(n)=mshift+1; %dummy them, if these criteria are met
%         else
    %         interpPGSSn=interp(sumsPGSSn(n,:),iup,3);
    %         interpPGSIn=interp(sumsPGSIn(n,:),iup,3);
    %         interpSISSn=interp(sumsSISSn(n,:),iup,3);
    %         leninterp=length(interpPGSSn);
    %         [xcmaxinterpPGSSn,imaxinterpPGSS]=max(interpPGSSn(1:leninterp-(iup-1)));
    %         [xcmaxinterpPGSIn,imaxinterpPGSI]=max(interpPGSIn(1:leninterp-(iup-1)));
    %         [xcmaxinterpSISSn,imaxinterpSISS]=max(interpSISSn(1:leninterp-(iup-1)));
    %         xcmaxconprev=-99999.;  %used to be 0; not good with glitches
    %         for iPGSS=max(1,imaxinterpPGSS-3*iup):min(imaxinterpPGSS+3*iup,iup*(2*mshift+1)-(iup-1)) %3 samples from peak; 
    %                                                                                  %intentionally wider than acceptable;
    %                                                                                  %iup-1 are extrapolated points
    %             for iPGSI=max(1,imaxinterpPGSI-3*iup):min(imaxinterpPGSI+3*iup,iup*(2*mshift+1)-(iup-1))
    %                 ibangon = (iup*mshift+1)-iPGSI+iPGSS;
    %                 if ibangon >= 1 && ibangon<=iup*(2*mshift+1)
    %                     xcmaxcon=interpPGSSn(iPGSS)+interpPGSIn(iPGSI)+interpSISSn(ibangon);
    %                     if xcmaxcon > xcmaxconprev
    %                         xcmaxconprev=xcmaxcon;
    %                         iPGSSbang=iPGSS;
    %                         iPGSIbang=iPGSI;
    %                     end
    %                 end
    %             end
    %         end
            xcmaxconprev=-99999.;  %used to be 0; not good with glitches
            imaxSTA12n=imaxSTA12(n); %This "n" for nth window; other "n's" for "normalized".  Unfortunately.
            imaxSTA13n=imaxSTA13(n);
            imaxSTA32n=imaxSTA32(n);
            sumsSTA12nn=sumsSTA12n(n,:);
            sumsSTA13nn=sumsSTA13n(n,:);
            sumsSTA32nn=sumsSTA32n(n,:);
            for iSTA12 =     max(1,imaxSTA12n-floor(loopoffmax+1)):min(imaxSTA12n+floor(loopoffmax+1),2*mshift+1)
                for iSTA13 = max(1,imaxSTA13n-floor(loopoffmax+1)):min(imaxSTA13n+floor(loopoffmax+1),2*mshift+1)
                    ibangon = (mshift+1)-iSTA13+iSTA12;
                    if ibangon >= 1 && ibangon <= 2*mshift+1
                        xcmaxcon=sumsSTA12nn(iSTA12)+sumsSTA13nn(iSTA13)+sumsSTA32nn(ibangon);
                        if xcmaxcon > xcmaxconprev
                            xcmaxconprev=xcmaxcon;
                            iSTA12bang=iSTA12;
                            iSTA13bang=iSTA13;
                        end
                    end
                end
            end
    %         iSISSbang=(iup*mshift+1)-iPGSIbang+iPGSSbang;
    %         if abs(iPGSSbang-imaxinterpPGSS) <= loopoffmax*iup && ...
    %            abs(iPGSIbang-imaxinterpPGSI) <= loopoffmax*iup && ...
    %            abs(iSISSbang-imaxinterpSISS) <= loopoffmax*iup && ...
    %            interpPGSSn(iPGSSbang)+interpPGSIn(iPGSIbang)+interpSISSn(iSISSbang) >= 3*xcmaxAVEnmin
    %             xmaxPGSSntmp(n)=(iPGSSbang-(iup*mshift+1))/iup;
    %             xmaxPGSIntmp(n)=(iPGSIbang-(iup*mshift+1))/iup;
    %             xmaxSISSntmp(n)=(iSISSbang-(iup*mshift+1))/iup;
            iSTA32bang=(mshift+1)-iSTA13bang+iSTA12bang;
%             if abs(iSTA12bang-imaxSTA12n) <= loopoffmax && ...  %not sure if these 3 lines are satisfied automatically ...
%                abs(iSTA13bang-imaxSTA13n) <= loopoffmax && ...
%                abs(iSTA32bang-imaxSTA32n) <= loopoffmax && ...
%                sumsSTA12n(n,iSTA12bang)+sumsSTA13n(n,iSTA13bang)+sumsSTA32n(n,iSTA32bang) >= 3*xcmaxAVEnmin
                xmaxSTA12ntmp(n)=iSTA12bang-(mshift+1); %without interpolation this is just centering.
                xmaxSTA13ntmp(n)=iSTA13bang-(mshift+1);
                xmaxSTA32ntmp(n)=iSTA32bang-(mshift+1);

                %for plotting traces
                imaxSTA12wr=round(xmaxSTA12ntmp(n)); %without interpolation this is not needed.
                imaxSTA13wr=round(xmaxSTA13ntmp(n));
    % 
                    STA12tr=STAopt(1,istart(n):iend(n)).*STAopt(2,istart(n)-imaxSTA12wr:iend(n)-imaxSTA12wr);
                    STA13tr=STAopt(1,istart(n):iend(n)).*STAopt(3,istart(n)-imaxSTA13wr:iend(n)-imaxSTA13wr);
                    STA32tr=STAopt(3,istart(n)-imaxSTA13wr:iend(n)-imaxSTA13wr).*STAopt(2,istart(n)-imaxSTA12wr:iend(n)-imaxSTA12wr);
                    cumsumtr=cumsum(STA12tr)+cumsum(STA13tr)+cumsum(STA32tr);
                    [cumsumtrdiff idiff]=max(cumsumtr(41:winlen)-cumsumtr(1:winlen-40));

                    STA1file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart(n):iend(n))' STAopt(1,istart(n):iend(n))'];
                    STA2file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart(n):iend(n))' STAopt(2,istart(n)-imaxSTA12wr:iend(n)-imaxSTA12wr)'];
                    STA3file(nin*winlen+1:(nin+1)*winlen,1:2)=[timsSTA(istart(n):iend(n))' STAopt(3,istart(n)-imaxSTA13wr:iend(n)-imaxSTA13wr)'];
                    STA12file(nin+1,1:2)=[imaxSTA12wr xcmaxSTA12n(n)];
                    STA13file(nin+1,1:2)=[imaxSTA13wr xcmaxSTA13n(n)];
                    STA32file(nin+1,1:3)=[cumsumtrdiff/cumsumtr(winlen) xcmaxSTA32n(n) idiff];

                    nin=nin+1;
                    istartkeep(nin)=istart(n); %For adding other stations later
                    iendkeep(nin)=iend(n); %For adding other stations later
                    aSTA12keep(nin,:)=[timswin(n) aSTA12(n)];
                    aSTA13keep(nin,:)=[timswin(n) aSTA13(n)];
                    aSTA32keep(nin,:)=[timswin(n) aSTA32(n)];
                    loopoffkeep(nin,:)=[timswin(n) loopoff(n)];
                    mapfile(nin,:)=[timswin(n) xmaxSTA13ntmp(n) xmaxSTA12ntmp(n) ...
                        xcmaxAVEn(n) loopoff(n) AmpComp(n) cumsumtrdiff timswin(n)-2+idiff/40. cumsumtrdiff/cumsumtr(winlen)];
%                else
%                    xmaxSTA12ntmp(n)=20; xmaxSTA13ntmp(n)=20; xmaxSTA32ntmp(n)=20;
%                end
%             else
%                 xmaxSTA12ntmp(n)=20; xmaxSTA13ntmp(n)=20; xmaxSTA32ntmp(n)=20; 
%             end
%         end
    end
    %fid = fopen(['ARMMAP/MAPS/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fid = fopen(['ARMMAP/map',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    fprintf(fid,'%9.1f %9.5f %9.5f %9.4f %9.2f %9.3f %10.4f %10.4f %8.3f\n',mapfile(1:nin,:)');
    fclose(fid);
    figure 
    subplot(4,1,1,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timsSTA(winlen:2*winlen),7+zeros(winlen+1,1),'k','linewidth',2);
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([0 timbig/2 -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,2,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([timbig/2 timbig -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,3,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([timbig 3*timbig/2 -mshift mshift]);
    ylabel('samples')
    box on
    subplot(4,1,4,'align'); 
    hold on
    plot(timswin,xcmaxAVEnmin*mshift+zeros(nwin,1),'k:');
    plot(timswin,zeros(nwin,1),'k:');
    plot(timswin,xcmaxAVEn*mshift,'g');
    plot(timswin,xmaxSTA12ntmp,'bs','MarkerSize',2);
    plot(timswin,xmaxSTA13ntmp,'ro','MarkerSize',2);
    plot(timswin,xmaxSTA32ntmp,'k*','MarkerSize',2);
    axis([3*timbig/2 2*timbig -mshift mshift]);
    xlabel('sec')
    ylabel('samples')
    box on
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    orient landscape
    % 
    figure
    colormap(jet)
    scatter(xmaxSTA13n-xmaxSTA12n+xmaxSTA32n,xcmaxAVEn,3,AmpComp)
    hold on 
    plot(-50:50,xcmaxAVEnmin+zeros(101,1),'k:');
    axis([min(-5,-2.5*loopoffmax) max(5,2.5*loopoffmax) -0.2 1.0])
    hrf = plotreflinesr(gca,-loopoffmax,'x','k');colorbar
    hrf = plotreflinesr(gca,loopoffmax,'x','k');colorbar
    box on
    title([IDENTIF,'_{',num2str(lo),'-',num2str(hi),'}'])
    print('-depsc',['FIGS/',IDENTIF,'-',int2str(winlen/40),'s_',num2str(lo),'-',num2str(hi),'e.eps'])

    
    % Add additional stations
    clear STAopt
    clear STAort
    clear STAnzeros
%     stasnew=['LZB '
    stasnew=['SSIB '
             'SILB '
             'SILBz'];
    nstanew=size(stasnew,1);
    STAoptnew=zeros(nstanew,nin*winlen);
    hi=12.;
    lo=1.5;
    for istanew=1:nstanew
        found=0;
        [LIA,idx]=ismember(stasnew(istanew,:),PERMSTA,'rows');
        if LIA
            found=found+LIA;
            if strcmp(PERMSTA(idx,1:3),'PGC')
                fact=1.6e-4;
            elseif strcmp(PERMSTA(idx,1:3),'LZB')
                fact=4.e-3;
            end
            [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        [LIA,idx]=ismember(stasnew(istanew,:),POLSTA,'rows');
        if LIA
            found=found+LIA; %better be 1
            if year==2003 && jday<213
                fact=20.0e-3;
            else
                fact=4.0e-3; 
            end
            [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
        end
        found=found
        for n=1:nin
%             istart=istartkeep(n); 
%             iend=istart+winlen-1;
            STAoptnew(istanew,(n-1)*winlen+1:n*winlen)=opt(istartkeep(n):iendkeep(n))';
        end
    end
    in=zeros(istanew,nin);
    ioff=zeros(istanew,nin);
    ccmax=zeros(istanew,nin);
    loff=zeros(istanew,nin);
    mshiftnew=ceil(1.6*mshift); %1.8
    for istanew=1:nstanew
        for n=1:nin
            is=(n-1)*winlen+1+iccwinoff;
            ie=is+winlenCC-1;
            cc1=xcorr(STA1file(is:ie,2),STAoptnew(istanew,is:ie),mshiftnew,'coeff');
            cc2=xcorr(STA2file(is:ie,2),STAoptnew(istanew,is:ie),mshiftnew,'coeff');
            cc3=xcorr(STA3file(is:ie,2),STAoptnew(istanew,is:ie),mshiftnew,'coeff');
            [cc1max,icc1]=max(cc1);
            [cc2max,icc2]=max(cc2);
            [cc3max,icc3]=max(cc3);
            ccmax(istanew,n)=(cc1max+cc2max+cc3max)/3;
            loff(istanew,n)=(abs(icc1-icc2)+abs(icc1-icc3)+abs(icc2-icc3))/3;
            if abs(icc1-icc2)<=20 && abs(icc1-icc3)<=20 && abs(icc2-icc3)<=20 && ...
                    (cc1max+cc2max+cc3max)/3>0.1*xcmaxAVEnmin %0.8
                in(istanew,n)=1;
                ioff(istanew,n)=round((cc1max*icc1+cc2max*icc2+cc3max*icc3)/(cc1max+cc2max+cc3max))-(1+mshiftnew); %-(1+mshiftnew) zeros it
                %Now shift new station (for plotting purposes) FOR HORIZONTALS ONLY
                if stasnew(istanew,5:5)==' '
                    if ioff(istanew,n)>0
                        STAoptnew(istanew,is+ioff(istanew,n):ie)=STAoptnew(istanew,is:ie-ioff(istanew,n));
                        STAoptnew(istanew,is:is+ioff(istanew,n)-1)=0;
                    elseif ioff(istanew,n)<0
                        STAoptnew(istanew,is:ie+ioff(istanew,n))=STAoptnew(istanew,is-ioff(istanew,n):ie);
                        STAoptnew(istanew,ie+ioff(istanew,n)+1:ie)=0;                    
                    end
                end
            end
        end
    end

    if winlen<=500
    scrsz=get(0,'ScreenSize');
    nt=0;
    nrow=4;
    mcol=3;
    for ifig=1:floor(nin/(nrow*mcol))+1
        figure('Position',[scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 9*scrsz(4)/10]);
        for n = 1:nrow
            for m = 1:mcol
                nt=nt+1;
                if nt <= nin
                     %if STA12file(nt,1) >= 10 && STA12file(nt,1) <= 16 && STA13file(nt,1) >= 2 && STA13file(nt,1) <= 8
                        subplot(3*nrow,mcol,3*(n-1)*mcol+m,'align');
                        plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STA1file(winlen*(nt-1)+1:winlen*nt,2),'r')
                        hold on
                        plot(STA2file(winlen*(nt-1)+1:winlen*nt,1),STA2file(winlen*(nt-1)+1:winlen*nt,2),'b')
                        plot(STA3file(winlen*(nt-1)+1:winlen*nt,1),STA3file(winlen*(nt-1)+1:winlen*nt,2),'k')
                        is = STA1file(winlen*(nt-1)+1,1);
                        ien= STA1file(winlen*nt,1);
                        %yma=0.4;
                        yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        xvect=[is is+2*(yma-ymi)*(winlen/160.)]; %amplitude bar originally scaled for 4-s window
                        yma=2.4*max(yma,-ymi);
                        yvect=[-0.9*yma -0.9*yma];
                        plot(xvect,yvect,'r','linewidth',3)
                        plot([is+1/hi is+1/lo],[-0.8*yma -0.8*yma],'k','linewidth',3)
                        text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
                        text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
                        istanew=1;
                        if stasnew(istanew,5:5)=='z' 
                            ymaz=-2.4*max(abs(STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt))); %If z component, flip sign and enlarge.
                        else
                            ymaz=yma;
                        end
                        if sum(in(:,nt))>=1
                            plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),(yma/ymaz)*STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt),'c')
                            text(0.5*(is+ien)-0.2, 0.66*yma, int2str(ioff(istanew,nt)),'fontsize',6);
                            text(0.5*(is+ien)-0.8, -0.8*yma, num2str(ccmax(istanew,nt),3),'fontsize',6);
                            text(0.5*(is+ien)+0.4, -0.8*yma, num2str(loff(istanew,nt),3),'fontsize',6);
                        end
                        box on
                        axis([is ien -yma yma])
                        hrf = plotreflinesr(gca,mapfile(nt,1),'x','k');
                        set(gca,'XTick',[is (is+ien)/2],'fontsize',6);

                        subplot(3*nrow,mcol,3*(n-1)*mcol+mcol+m,'align');
                        plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STA1file(winlen*(nt-1)+1:winlen*nt,2),'r')
                        hold on
                        plot(STA2file(winlen*(nt-1)+1:winlen*nt,1),STA2file(winlen*(nt-1)+1:winlen*nt,2),'b')
                        plot(STA3file(winlen*(nt-1)+1:winlen*nt,1),STA3file(winlen*(nt-1)+1:winlen*nt,2),'k')
                        is = STA1file(winlen*(nt-1)+1,1);
                        ien= STA1file(winlen*nt,1);
                        %yma=0.4;
                        yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        yma=2.4*max(yma,-ymi);
                        text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
                        text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
                        istanew=2;
                        if stasnew(istanew,5:5)=='z' 
                            ymaz=-2.4*max(abs(STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt))); %If z component, flip sign and enlarge.
                        else
                            ymaz=yma;
                        end
                        if sum(in(:,nt))>=1
                            plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),(yma/ymaz)*STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt),'c')
                            text(0.5*(is+ien)-0.2, 0.66*yma, int2str(ioff(istanew,nt)),'fontsize',6);
                            text(0.5*(is+ien)-0.8, -0.8*yma, num2str(ccmax(istanew,nt),3),'fontsize',6);
                            text(0.5*(is+ien)+0.4, -0.8*yma, num2str(loff(istanew,nt),3),'fontsize',6);
                        end
                        box on
                        axis([is ien -yma yma])
                        set(gca,'XTick',[is (is+ien)/2],'fontsize',6);

                        subplot(3*nrow,mcol,3*(n-1)*mcol+2*mcol+m,'align');
                        plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),STA1file(winlen*(nt-1)+1:winlen*nt,2),'r')
                        hold on
                        plot(STA2file(winlen*(nt-1)+1:winlen*nt,1),STA2file(winlen*(nt-1)+1:winlen*nt,2),'b')
                        plot(STA3file(winlen*(nt-1)+1:winlen*nt,1),STA3file(winlen*(nt-1)+1:winlen*nt,2),'k')
                        is = STA1file(winlen*(nt-1)+1,1);
                        ien= STA1file(winlen*nt,1);
                        yma=max(max([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        ymi=min(min([STA1file(winlen*(nt-1)+1:winlen*nt,2) STA2file(winlen*(nt-1)+1:winlen*nt,2) ...
                            STA3file(winlen*(nt-1)+1:winlen*nt,2)]));
                        yma=2.4*max(yma,-ymi);
                        text(is+0.2, 0.66*yma, int2str(STA13file(nt,1)),'fontsize',6);
                        text(ien-0.6, 0.66*yma, int2str(STA12file(nt,1)),'fontsize',6);
                        istanew=3;
                        if stasnew(istanew,5:5)=='z' 
                            ymaz=-2.4*max(abs(STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt))); %If z component, flip sign and enlarge.
                        else
                            ymaz=yma;
                        end
                        if sum(in(:,nt))>=1
                            plot(STA1file(winlen*(nt-1)+1:winlen*nt,1),(yma/ymaz)*STAoptnew(istanew,winlen*(nt-1)+1:winlen*nt),'c')
                            text(0.5*(is+ien)-0.2, 0.66*yma, int2str(ioff(istanew,nt)),'fontsize',6);
                            text(0.5*(is+ien)-0.8, -0.8*yma, num2str(ccmax(istanew,nt),3),'fontsize',6);
                            text(0.5*(is+ien)+0.4, -0.8*yma, num2str(loff(istanew,nt),3),'fontsize',6);
                        end
                        box on
                        axis([is ien -yma yma])
                        hrf = plotreflinesr(gca,mapfile(nt,1)-5.064,'x','k');
                        set(gca,'XTick',[is (is+ien)/2],'fontsize',6);
                end
            end
        end
        orient landscape
    %     if ifig <= 9
    %         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'.',int2str(0),int2str(0),int2str(ifig),'.eps'])
    %     elseif ifig <= 99
    %         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'.',int2str(0),int2str(ifig),'.eps'])
    %     else
    %         print('-depsc',['ARMMAP/WIGS/',IDENTIF,'-',num2str(lo),'-',num2str(hi),'.',int2str(ifig),'.eps'])
    %     end

    end
    % %fid = fopen(['ARMMAP/MAPS/pks',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    % fid = fopen(['ARMMAP/pks',IDENTIF,'_',num2str(lo),'-',num2str(hi),'-','ms',int2str(mshift),'-',int2str(winlen/40),'s'],'w');
    % fprintf(fid,'%9.1f %9.3f %10.6f %9.3f %9.3f %9.3f %9.3f \n',pkfile(1:nin,:)');
    % fclose(fid);
    end

medlok=median(abs(loopoffkeep))
medaSTA12=median(aSTA12keep)
medaSTA13=median(aSTA13keep)
medaSTA32=median(aSTA32keep)

end


