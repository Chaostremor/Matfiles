% travel time curve fitting with Pb V9
% 
% USAGE:
%     construct a four-layered model, sediment, upper crust, lower crust and mantle,
%     which the source is located in the upper crust
%     use the known and unknown parameters to calculate the theoratical curve and fit the real curve best
% Model:     (12)
%      tk1, vp1, vs1
%      source
%      tk2, vp2, vs2
%      tk3, vp3, vs3
%            vp4, vs4
%      
%
%  Changes against last version:
%      Consider the pws result of Pb is solid, the lower crust velocity is
%      known too, compared to v7, use stack pn arrival and fit pg arrival, some kind of average 
%
% Author: C. Song,  2017.5.27

clear; 

%% 1. settle up the models
% load manually picked pg and sg arrival times named 'pg' and 'sg' 
load('pgsgpicktime.mat');
ndist = length(dist);           % num. of data

% results from PWS 600-400 result in V6, refdist=498.6431
% vpn = 8.40;                          % vel. of Pn
% tp0 = 72.12;                         % time of Pn at refdist 498.6431
% dtp1 = 5.24;                         % difference time between Pn and pPn
% dtp2 = 7.57;                         % difference time between Pn and sPn
% vsn = 4.61;                           % vel. of Sn
% ts0 = 122.8;                          % time of Sn at refdist 498.6431
% dts1 = 8.0;                            % difference time between Sn and sSn
% refdist = 498.6431;

% results from PWS 400-250 result in V6, 
vpn = 8.12;
perip = 1.66;                        % 1/4 period length of stacked pn phase , 50.54-48.88   
tp0 = 50.54-perip;                % time of Pn at refdist - 1/4 period, to make the time from peak to the first motion
dtp1 = 5.94;                         % difference time between Pn and pPn
% dtp2 = 7.57;                         % difference time between Pn and sPn
vpb = 6.75;                          % vel. of Pb
peripb = 1.02;                      % 1/4 period length of stacked pn phase , 57.75-56.73  
tp1 = 57.75-peripb;              % time of Pb at refdist - 1/4 period, to make the time from peak to the first motion
vsn = 4.48;                           % vel. of Sn
peris = 1.46;                         % 1/4 period length of stacked sn phase , 85.96-84.5
ts0 = 85.96-peris;                  % time of Sn at refdist - 1/4 period, to make the time from peak to the first motion
dts1 = 8.0;                            % difference time between Sn and sSn
refdist = 326.0027;

% known parameters (4)
%vp2 = 5.87;  vp2 = 5.467
%vs2 = 3.48;
%vp4 = vpn;
%vs4 = vsn;

a1 = 6.188;
b1 = 33.83;
distcurve = 0:0.1:400;                                         % continuous sampling of dist
pgcurve = a1.*sqrt((distcurve./b1).^2+1);           % pg travel time picks fitting curve
pncurve = (distcurve - refdist)/vpn + tp0;          % PWS pn curve
pnstack = (dist- refdist)/vpn + tp0;                            % dist is real data,
pbcurve = (distcurve - refdist)/vpb + tp1;    % PWS pn curve
pbstack = (dist- refdist)/vpb + tp1;


a2 = 10.01;
b2 = 34.37;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);            % sg travel time picks fitting curve
sncurve = (distcurve - refdist)/vsn + ts0;             % PWS sn curve   
snstack = (dist - refdist)/vsn + ts0;

% unknown ones  (7)
vp1 = 3.0: 0.1: 5.0;
nvp1 = length(vp1);

tk1 = 3.5: 0.1: 6;
ntk1 = length(tk1);           % vp1, tk1:     sediment layer        

% vp2 = 5.4: 0.1: 5.8;
% vp2 = 5.65: 0.1: 5.85;
% vp2 = 5.5: 0.02: 5.9;
vp2 = 5.0: 0.02: 5.6;
nvp2 = length(vp2);

tk2 = 11: 0.2: 19;
ntk2 = length(tk2);           % vp2, tk2:     upper source layer of upper crust 

tk3 = 11: 0.2: 19;
ntk3 = length(tk3);           % vp2, tk3:     lower source layer of upper crust

% vp4 = 6.1: 0.2: 8.3;
% nvp4 = length(vp4);

tk4 = 10: 0.4: 25;
ntk4 = length(tk4);          % tk4:     lower crust

% vs1 = 1.0: 0.2: 3.4;
% nvs1 = length(vs1);
% 
% vs3 = 3.5: 0.2: 4.6;
% nvs3 = length(vs3);


%% 2. get the speed for P wave and thickness of each layer (7)
%raypg = 0.1680: 0.0004: 0.1696;       % 1/5.87 ~= 0.1703,  max(raypg) < 1/max(vp2)
% raypg = 0.1688: 0.0001: 0.1692;
raypg = 0.16235: 0.0002: 0.18235;
nraypg = length(raypg);
raypn = 1/vpn;
raypb = 1/vpb;

% 2.1 sediment layer
% pn
size1 = nvp1*ntk1;  
div1pn = zeros(size1, 1);
mul1pn = zeros(size1, 1);
ict = 1;
for ii = 1:nvp1
    etapn = sqrt((1.0/vp1(ii))^2-raypn^2);
    for jj = 1:ntk1            
        div1pn(ict) = tk1(jj)/etapn;
        mul1pn(ict) = tk1(jj)*etapn;        
        ict = ict+1;
    end
end
% pg
div1pg = zeros(size1, nraypg);
mul1pg = zeros(size1, nraypg);
for mm = 1:nraypg
    ict = 1;
    for ii = 1:nvp1
        etapg = sqrt((1.0/vp1(ii))^2-raypg(mm)^2);
        for jj = 1:ntk1            
            div1pg(ict, mm) = tk1(jj)/etapg;
            mul1pg(ict, mm) = tk1(jj)*etapg;        
            ict = ict+1;
        end
    end
end
% pb
div1pb = zeros(size1, 1);
mul1pb = zeros(size1, 1);
ict = 1;
for ii = 1:nvp1
    etapb = sqrt((1.0/vp1(ii))^2-raypb^2);
    for jj = 1:ntk1            
        div1pb(ict) = tk1(jj)/etapb;
        mul1pb(ict) = tk1(jj)*etapb;        
        ict = ict+1;
    end
end

% 2.2 upper crust 
% pn in upper source layer 
size2 = nvp2*ntk2; 
div2pn = zeros(size2, 1);
mul2pn = zeros(size2, 1);     
% pn in lower source layer
size3 = nvp2*ntk3;
div3pn = zeros(size3, 1);
mul3pn = zeros(size3, 1);   
ict1 = 1;
ict2 = 1;
for ii = 1:nvp2
    etapn = sqrt((1.0/vp2(ii))^2-raypn^2);
    % upper source layer 
    for jj = 1:ntk2        
        div2pn(ict1) = tk2(jj)/etapn;
        mul2pn(ict1) = tk2(jj)*etapn;
        ict1 = ict1+1;
    end
    % lower source layer 
    for jj = 1:ntk3
        div3pn(ict2) = tk3(jj)/etapn;
        mul3pn(ict2) = tk3(jj)*etapn;
        ict2 = ict2+1;
    end
end
% pg in upper source layer 
div2pg = zeros(size2, nraypg);
mul2pg = zeros(size2, nraypg);
for mm = 1:nraypg
    ict1 = 1;
    for ii = 1:nvp2
        etapg = sqrt((1.0/vp2(ii))^2-raypg(mm)^2);
        for jj = 1:ntk2        
            div2pg(ict1, mm) = tk2(jj)/etapg;
            mul2pg(ict1, mm) = tk2(jj)*etapg;
            ict1 = ict1+1;
        end
    end
end
% pb in upper source layer 
div2pb = zeros(size2, 1);
mul2pb = zeros(size2, 1);     
% pb in lower source layer
div3pb = zeros(size3, 1);
mul3pb = zeros(size3, 1);  
ict1 = 1;
ict2 = 1;
for ii = 1:nvp2
    etapb = sqrt((1.0/vp2(ii))^2-raypb^2);
    % upper source layer 
    for jj = 1:ntk2        
        div2pb(ict1) = tk2(jj)/etapb;
        mul2pb(ict1) = tk2(jj)*etapb;
        ict1 = ict1+1;
    end
    % lower source layer 
    for jj = 1:ntk3
        div3pb(ict2) = tk3(jj)/etapb;
        mul3pb(ict2) = tk3(jj)*etapb;
        ict2 = ict2+1;
    end
end

% 2.3 lower crust
% pn
div4pn = zeros(ntk4, 1);
mul4pn = zeros(ntk4, 1);      
etapn = sqrt((1.0/vpb)^2-raypn^2);
for jj = 1:ntk4         
    div4pn(jj) = tk4(jj)/etapn;
    mul4pn(jj) = tk4(jj)*etapn; 
end

maxdist = dist(end);
mindist = dist(1);
minmisp = 10000;
minii = 1;
minjj = 1;
minkk = 1;
minmm = 1;
tpn0pws = tp0 - raypn*refdist;
tpb0pws = tcd p1 - raypb*refdist;
for ii = 1: size1
    ii
    for jj = 1: size2
        
        % part 1, Pg(mod) ~= Pg(fit)
        distuse = [];
        pgmuse = [];
        pgduse = [];
        [~, ivp2] = ind2sub([ntk2, nvp2], jj);
        for nn = 1: nraypg
            if (1/vp2(ivp2) < raypg(nn))
                continue;
            end
            x = raypg(nn)*(div1pg(ii, nn)+div2pg(jj, nn));
            if (x < 0 || x > maxdist)
                continue;
            end
            pgmod = raypg(nn)*x+mul1pg(ii, nn)+mul2pg(jj, nn);
            pgfit = a1.*sqrt((x./b1).^2+1);
            distuse = [distuse; x];
            pgmuse = [pgmuse; pgmod];
            pgduse = [pgduse; pgfit];            
        end
        if isempty(distuse)
            continue;
        end
        npg = length(distuse);
        mis1 = sum((pgmuse-pgduse).^2)/npg;
                
        % part 2, pPn-Pn(mod) ~= pPn-Pn(stack)
        dtpmod = 2*(mul1pn(ii)+mul2pn(jj));
        mis2 = (dtpmod- dtp1).^2;
        
        for kk = 1: ntk3
            
            % part 5, Pb intercept with x = 0
            indice = (ivp2-1)*ntk3 + kk;  
            xmin = raypb*(div1pb(ii)+div2pb(jj)+2*div3pb(indice));
            tmin = raypb*xmin + mul1pb(ii)+mul2pb(jj)+2*mul3pb(indice);
            tpb0mod = mul1pb(ii)+mul2pb(jj)+2*mul3pb(indice);
            mis5 = (tpb0mod-tpb0pws).^2;
                
            % part 6, Pb-Pg
            pbmuse = raypb*(distuse-xmin)+tmin;
            dtmod = pbmuse - pgmuse;
            pbduse = (distuse- refdist)/vpb + tp1;
            dtdata = pbduse -pgduse;
            mis6 = sum((dtmod - dtdata).^2)/npg;
                
            for mm = 1: ntk4
                                
                % part 3, Pn intercept with x = 0    
                xmin = raypn*(div1pn(ii)+div2pn(jj)+2*div3pn(indice)+2*div4pn(mm));
                tmin = raypn*xmin + mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                tpn0mod = mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                mis3 = (tpn0mod-tpn0pws).^2;
                
                % part 4, Pn-Pg
                pnmuse = raypn*(distuse-xmin)+tmin;
                dtmod = pnmuse - pgmuse;
                pnduse = (distuse- refdist)/vpn + tp0;
                dtdata = pnduse -pgduse;
                mis4 =sum((dtmod - dtdata).^2)/npg;
               
                
                % sum all parts
                misfit = mis1+mis2+mis3+mis4+mis5+mis6;
                if misfit < minmisp
                    minmisp = misfit;
                    minmis1 = mis1;
                    minmis2 = mis2;
                    minmis3 = mis3;
                    minmis4 = mis4;
                    minmis5 = mis5;
                    minmis6 = mis6;
                    minii = ii;
                    minjj = jj;
                    minkk = kk;
                    minmm = mm;
                end
            end
        end                
    end
end
        
[itk1, ivp1] = ind2sub([ntk1, nvp1], minii);
[itk2, ivp2] = ind2sub([ntk2, nvp2], minjj);
itk3 = minkk;
itk4 = minmm;
bestvp1 = vp1(ivp1);
bestvp2 = vp2(ivp2);
besttk1 = tk1(itk1);
besttk2 = tk2(itk2);
besttk3 = tk3(itk3);
besttk4 = tk4(itk4);
bestdep = besttk1+besttk2;        
        
save('p_wave_para.mat');      


%% 3. get the speed for S wave of each layer (3)
% set para. range 
vs1 = 1.5: 0.1: 2.9;
nvs1 = length(vs1);

vs2 = 3.3: 0.02: 3.7;
nvs2 = length(vs2);

vs4 = 3.5: 0.1: 4.3;
nvs4 = length(vs4);

raysn =  1/vsn;
raysg = 0.27028: 0.0002: 0.29028;    % P(x)=[0, 0.2902], 
nraysg = length(raysg);

% 3.1 sediment layer
% sn
div1sn = zeros(nvs1, 1);
mul1sn = zeros(nvs1, 1);      
for ii = 1:nvs1
    etasn = sqrt((1.0/vs1(ii))^2-raysn^2);
    div1sn(ii) = besttk1/etasn;
    mul1sn(ii) = besttk1*etasn;      
end
% sg
div1sg = zeros(nvs1, nraysg);
mul1sg = zeros(nvs1, nraysg);
for mm = 1: nraysg
    for ii = 1:nvs1
        etasg = sqrt((1.0/vs1(ii))^2-raysg(mm)^2);
        div1sg(ii, mm) = besttk1/etasg;
        mul1sg(ii, mm) = besttk1*etasg;  
    end
end

% 3.2 upper crust
% sn in upper source layer
div2sn = zeros(nvs2, 1);
mul2sn = zeros(nvs2, 1);      
% sn in lower source layer
div3sn = zeros(nvs2, 1);
mul3sn = zeros(nvs2, 1);    
for ii = 1:nvs2
    etasn = sqrt((1.0/vs2(ii))^2-raysn^2);
    % upper source layer
    div2sn(ii) = besttk2/etasn;
    mul2sn(ii) = besttk2*etasn;
    % lower source layer
    div3sn(ii) = besttk3/etasn;
    mul3sn(ii) = besttk3*etasn; 
end
% sg in upper source layer
div2sg = zeros(nvs2, nraysg);
mul2sg = zeros(nvs2, nraysg);
for mm = 1: nraysg
    for ii = 1:nvs2
        etasg = sqrt((1.0/vs2(ii))^2-raysg(mm)^2);
        div2sg(ii, mm) = besttk2/etasg;
        mul2sg(ii, mm) = besttk2*etasg;  
    end
end

% 3.3 lower crust
% sn
div4sn = zeros(nvs4, 1);
mul4sn = zeros(nvs4, 1);      
for ii = 1:nvs4
    etasn = sqrt((1.0/vs4(ii))^2-raysn^2);
    div4sn(ii) = besttk4/etasn;
    mul4sn(ii) = besttk4*etasn; 
end

minmiss = 10000;
minii = 1;
minjj = 1;
minkk = 1;
tsn0pws = ts0 - raysn*refdist;
for ii = 1: nvs1
    ii
    for jj = 1: nvs2
        
        % part 1, Sg(mod) ~= Sg(pick)
        distuse = [];
        sgmuse = [];
        sgduse = [];
       for nn = 1: nraysg
            if (1/vs2(jj) < raysg(nn))
                continue;
            end
            x = raysg(nn)*(div1sg(ii, nn)+div2sg(jj, nn));
            if (x < 0 || x > maxdist)
                continue;
            end
            sgmod = raysg(nn)*x+mul1sg(ii, nn)+mul2sg(jj, nn);
            sgfit = a2.*sqrt((x./b2).^2+1);
            distuse = [distuse; x];
            sgmuse = [sgmuse; sgmod];
            sgduse = [sgduse; sgfit];
        end
        if isempty(distuse)
            continue;
        end
        nsg = length(distuse);
        mis1 =sum((sgmuse-sgduse).^2)/nsg;

        for kk = 1: nvs4

            % part 2, Sn intercept with x = 0
            xmin = raysn*(div1sn(ii)+div2sn(jj)+2*div3sn(jj)+2*div4sn(kk));
            tmin = raysn*xmin + mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            tsn0mod = mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            mis2 = (tsn0mod-tsn0pws).^2;

            % part 3, Sn(stack) - Sg(fit)
            snmuse = raysn*(distuse-xmin)+tmin;
            dtmod = snmuse - sgmuse;
            snduse = (distuse- refdist)/vsn + ts0;
            dtdata = snduse -sgduse;
            mis3 = sum((dtmod - dtdata).^2)/nsg;

            % sum all parts
            misfit = mis1+mis2+mis3;
            if misfit < minmiss
                minmiss = misfit;
                minmis1 = mis1;
                minmis2 = mis2;
                minmis3 = mis3;
                minii = ii;
                minjj = jj;
                minkk = kk;
            end
        end                
    end
end
        
bestvs1 = vs1(minii);
bestvs2 = vs2(minjj);
bestvs4 = vs4(minkk);     
        
save('s_wave_para.mat');  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        