% travel time curve fitting V5
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
%      V4 also has some lethal flaws
%
% Author: C. Song,  2017.5.19

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
vpn = 8.1;
tp0 = 50.54;                         % time of Pn at refdist 
dtp1 = 5.94;                         % difference time between Pn and pPn
% dtp2 = 7.57;                         % difference time between Pn and sPn
vsn = 4.4;                           % vel. of Sn
ts0 = 86.15;                          % time of Sn at refdist 498.6431
dts1 = 8.0;                            % difference time between Sn and sSn
refdist = 326.0027;

% known parameters (4)
%vp2 = 5.87;  vp2 = 5.467
%vs2 = 3.48;   vs2 = 3.434
%vp4 = vpn;
%vs4 = vsn;

a1 = 6.188;
b1 = 33.83;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);           % pg travel time picks fitting curve
pncurve = (distcurve - refdist)/vpn + tp0;    % PWS pn curve
pn = (dist- refdist)/vpn + tp0;

a2 = 10.01;
b2 = 34.37;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
sncurve = (distcurve - refdist)/vsn + ts0;
sn = (dist - refdist)/vsn + ts0;

% unknown ones  (7)
vp1 = 2.5: 0.2: 5;                    
nvp1 = length(vp1);

tk1 = 0.4: 0.2: 6;
ntk1 = length(tk1);           % vp1, tk1:     sediment layer        

% vp2 = 5.4: 0.1: 5.8;
% vp2 = 5.65: 0.1: 5.85;
% vp2 = 5.5: 0.02: 5.9;
vp2 = 5.0: 0.02: 5.48;
nvp2 = length(vp2);      

tk2 = 5: 0.4: 20;
ntk2 = length(tk2);           % vp2, tk2:     upper source layer of upper crust 

tk3 = 5: 0.4: 15;             
ntk3 = length(tk3);           % vp2, tk3:     lower source layer of upper crust

vp4 = 6.1: 0.2: 8.0;
nvp4 = length(vp4);

tk4 = 10: 0.4: 25;
ntk4 = length(tk4);          % vp4, tk4:     lower crust

% vs1 = 1.0: 0.2: 3.4;
% nvs1 = length(vs1);
% 
% vs3 = 3.5: 0.2: 4.6;
% nvs3 = length(vs3);
      

%% 3. get the speed for S wave of each layer (3)
% set para. range 
vs1 = 1.5: 0.1: 3.0;
nvs1 = length(vs1);

vs2 = 3.0: 0.02: 3.44;
nvs2 = length(vs2);

vs4 = 3.5: 0.1: 4.3;
nvs4 = length(vs4);

raysn =  1/vsn;
raysg = 0.28425: 0.0015: 0.29025;    % P(x)=[0, 0.2902], 
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

minmis = 10000;
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
        temp = 0;
        nsg =zeros(nraysg, 1);
        for nn = 1: nraysg
            x = raysg(nn)*(div1sg(ii, nn)+div2sg(jj, nn));
            if (x < 100 || x > 400)
                nsg(nn) = 0;
                continue;
            end
            inddist = find((dist>=x-10) & (dist<=x+10));
            
%             inddist = find((dist>=x-10) & (dist<=x+10));
%             if isempty(inddist)
%                 npg(nn)=0;
%                 continue;
%             end
            t = raysg(nn)*x+mul1sg(ii, nn)+mul2sg(jj, nn);
            nsg(nn) = length(inddist);
            sgmod = raysg(nn)*(dist(inddist)-x)+t;
            temp = temp+sum((sgmod-sg(inddist)).^2)/nsg(nn);
            distuse = [distuse; dist(inddist)];
            sgmuse = [sgmuse; sgmod];
            sgduse = [sgduse; sg(inddist)];            
        end
        if isempty(distuse)
            continue;
        end
        indn = find(nsg == 0);
        nsg(indn) = [];
        snsg = length(nsg);
        mis1 = temp/snsg;
                        
        for kk = 1: nvs4
                                
            % part 2, Sn intercept with x = 0
            xmin = raysn*(div1sn(ii)+div2sn(jj)+2*div3sn(jj)+2*div4sn(kk));
            tmin = raysn*xmin + mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            tsn0mod = mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            mis2 = (tsn0mod-tsn0pws).^2;
                
            % part 3, Sn-Sg
            snmuse = raysn*(distuse-xmin)+tmin;
            dtmod = snmuse - sgmuse;
            snduse = (distuse- refdist)/vsn + ts0;
            dtdata = snduse -sgduse;
            dmd = dtmod - dtdata;
            temp = 0;
            for nn = 1: snsg
                inde = sum(nsg(1:nn));
                inds = inde-nsg(nn)+1;
                temp = temp+sum((dmd(inds: inde)).^2)/nsg(nn);
            end
            mis3 = temp/snsg;
                
            % sum all parts
            misfit = mis1+mis2+mis3;
            if misfit < minmis
                minmis = misfit;
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        