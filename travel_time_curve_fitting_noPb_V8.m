% travel time curve fitting without pb V8
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
%      use real Pn (from stack ) and Pg arrival (from pick) or Sn and Sg arrival
%
% Author: C. Song,  2017.5.26

clear; 

%% 1. settle up the models
% load manually picked pg and sg arrival times named 'pg' and 'sg' and real 'pn' and 'sn' arrival from free shift stack 
load('pgsgpicktime.mat');
ndist = length(dist);           % num. of data

% results from PWS 400-250 result in V6, 
vpn = 8.12;         
perip = 1.66;                        % 1/4 period length of stacked pn phase , 50.54-48.88   
tp0 = 50.54-perip;                % time of Pn at refdist - 1/4 period, to make the time from peak to the first motion
dtp1 = 5.94;                         % difference time between Pn and pPn
% dtp2 = 7.57;                         % difference time between Pn and sPn
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
pn = pn-perip;

a2 = 10.01;
b2 = 34.37;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);            % sg travel time picks fitting curve
sncurve = (distcurve - refdist)/vsn + ts0;             % PWS sn curve   
snstack = (dist - refdist)/vsn + ts0;
sn = sn-peris;

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

vp4 = 6.0: 0.1: 7.0;
nvp4 = length(vp4);

tk4 = 15: 0.2: 25;
ntk4 = length(tk4);          % vp4, tk4:     lower crust

% vs1 = 1.0: 0.2: 3.4;
% nvs1 = length(vs1);
% 
% vs3 = 3.5: 0.2: 4.6;
% nvs3 = length(vs3);


%% 2. get the speed for P wave and thickness of each layer (7)
%raypg = 0.1680: 0.0004: 0.1696;       % 1/5.87 ~= 0.1703,  max(raypg) < 1/max(vp2)
% raypg = 0.1688: 0.0001: 0.1692;
raypg = 0.17835: 0.001: 0.18235;
nraypg = length(raypg);
raypn = 1/vpn;

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

% 2.3 lower crust
% pn
size4 = nvp4*ntk4;
div4pn = zeros(size4, 1);
mul4pn = zeros(size4, 1);      
ict = 1;
for ii = 1:nvp4
    etapn = sqrt((1.0/vp4(ii))^2-raypn^2);
    for jj = 1:ntk4         
        div4pn(ict) = tk4(jj)/etapn;
        mul4pn(ict) = tk4(jj)*etapn; 
        ict = ict+1;
    end
end

maxdist = dist(end);
mindist = dist(1);
cutdist = dist(index);
minmisp = 10000;
minii = 1;
minjj = 1;
minkk = 1;
minmm = 1;
tpn0pws = tp0 - raypn*refdist;
for ii = 1: size1
    ii
    for jj = 1: size2
        
        % part 1, Pg(mod) ~= Pg(pick)
        distuse = [];
        pgmuse = [];
        pgduse = [];
        pnduse = [];
        temp = 0;
        npg = zeros(nraypg, 1);
        [~, ivp2] = ind2sub([ntk2, nvp2], jj);
        for nn = 1: nraypg
            if (1/vp2(ivp2) <= raypg(nn))                
              npg(nn)=0;
              continue;
            end
            x = raypg(nn)*(div1pg(ii, nn)+div2pg(jj, nn));
            if (x < cutdist || x > maxdist)                     % Pn has data only > 250 km ( >= cutdist)
                npg(nn) = 0;
                continue;
            end
            inddist = find((dist>=x-10) & (dist<=x+10) & (dist>=cutdist) & (dist<=maxdist));             % find index of dist which belongs [x-10, x+10]
            t = raypg(nn)*x+mul1pg(ii, nn)+mul2pg(jj, nn);         % cal t at x
            npg(nn) = length(inddist);                                         % cal num. of pg used 
            pgmod = raypg(nn)*(dist(inddist)-x)+t;                      % cal model pg time, dist and pg use same index
            temp = temp+sum((pgmod-pg(inddist)).^2)/npg(nn);
            
            % prepare for pn-pg
            distuse = [distuse; dist(inddist)];
            pgmuse = [pgmuse; pgmod];
            pgduse = [pgduse; pg(inddist)];
            pnduse = [pnduse; pn(inddist-index+1)];
        end
        if isempty(distuse)
            continue;
        end
        indn = find(npg == 0);
        npg(indn) = [];
        snpg = length(npg);
        mis1 = temp/snpg;
                
        % part 2, pPn-Pn(mod) ~= pPn-Pn(stack)
        dtpmod = 2*(mul1pn(ii)+mul2pn(jj));
        mis2 = (dtpmod- dtp1).^2;
        
        for kk = 1: ntk3
            for mm = 1: size4
                                
                % part 3, Pn(mod) intercept with x = 0 
                indice = (ivp2-1)*ntk3 + kk;      
                xmin = raypn*(div1pn(ii)+div2pn(jj)+2*div3pn(indice)+2*div4pn(mm));
                tmin = raypn*xmin + mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                tpn0mod = mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                mis3 = (tpn0mod-tpn0pws).^2;
                
                % part 4, Pn(data)-Pg(pick),  
                pnmuse = raypn*(distuse-xmin)+tmin;
                dtmod = pnmuse - pgmuse;
                dtdata = pnduse -pgduse;
                dmd = dtmod - dtdata;
                temp = 0;
                for nn = 1: snpg
                    inde = sum(npg(1: nn));
                    inds = inde-npg(nn)+1;
                    temp = temp+sum((dmd(inds: inde)).^2)/npg(nn);
                end
                mis4 = temp/snpg;
                
                % sum all parts
                misfit = mis1+mis2+mis3+mis4;
                if misfit < minmisp
                    minmisp = misfit;
                    minmis1 = mis1;
                    minmis2 = mis2;
                    minmis3 = mis3;
                    minmis4 = mis4;
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
[itk4, ivp4] = ind2sub([ntk4, nvp4], minmm);
bestvp1 = vp1(ivp1);
bestvp2 = vp2(ivp2);
bestvp4 = vp4(ivp4);
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
raysg = 0.28628: 0.001: 0.29028;    % P(x)=[0, 0.2902], 
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
        snduse = [];
        temp = 0;
        nsg =zeros(nraysg, 1);
        for nn = 1: nraysg
            if (1/vs2(jj) < raysg(nn))
                nsg(nn)=0;
                continue;
            end
            x = raysg(nn)*(div1sg(ii, nn)+div2sg(jj, nn));
            if (x < cutdist || x > maxdist)
                nsg(nn) = 0;
                continue;
            end
            inddist = find((dist>=x-10) & (dist<=x+10)& (dist>=cutdist) & (dist<=maxdist));
            t = raysg(nn)*x+mul1sg(ii, nn)+mul2sg(jj, nn);
            nsg(nn) = length(inddist);
            sgmod = raysg(nn)*(dist(inddist)-x)+t;
            temp = temp+sum((sgmod-sg(inddist)).^2)/nsg(nn);
            
            % prepare for sn-sg
            distuse = [distuse; dist(inddist)];
            sgmuse = [sgmuse; sgmod];
            sgduse = [sgduse; sg(inddist)];
            snduse = [snduse; sn(inddist+1-index)];
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