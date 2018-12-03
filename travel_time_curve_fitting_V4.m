% travel time curve fitting V4
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
%      V3 has some lethal flaws
%
% Author: C.Song,  2017.5.13

clear; 
%% 1. settle up the models
% load manually picked pg and sg arrival times named 'pg' and 'sg' 
load('pgsgpicktime.mat');
ndist = length(dist);           % num. of data

% results from PWS
vpn = 8.40;                          % vel. of Pn
tp0 = 72.12;                         % time of Pn at refdist 498.6431
dtp1 = 5.24;                         % difference time between Pn and pPn
dtp2 = 7.57;                         % difference time between Pn and sPn
vsn = 4.61;                           % vel. of Sn
ts0 = 122.8;                          % time of Sn at refdist 498.6431
dts1 = 8.0;                            % difference time between Sn and sSn

% known parameters (4)
%vp2 = 5.87;
%vs2 = 3.48;
%vp4 = vpn;
%vs4 = vsn;

a1 = 5.969;
b1 = 35.05;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);           % pg travel time picks fitting curve
pncurve = (distcurve - 498.6431)/vpn + tp0;    % PWS pn curve
pn = (dist- 498.6431)/vpn + tp0;

a2 = 10.34;
b2 = 35.95;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
sncurve = (distcurve - 498.6431)/vsn + ts0;
sn = (dist - 498.6431)/vsn + ts0;

% figure
% plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;
% plot(sg, dist, 'k.', 'MarkerSize', 8); hold on;
% plot(pgcurve, distcurve, 'b-'); hold on;
% plot(sgcurve, distcurve, 'b-'); hold on;
% plot(pncurve, distcurve, 'r-'); hold on;
% plot(sncurve, distcurve, 'r-'); hold on;

% unknown ones  (7)
vp1 = 2.5: 0.2: 5;                    
nvp1 = length(vp1);

tk1 = 0.1: 0.2: 5;
ntk1 = length(tk1);           % vp1, tk1:     sediment layer        

% vp2 = 5.4: 0.1: 5.8;
vp2 = 5.65: 0.1: 5.85;
nvp2 = length(vp2);      

tk2 = 5: 0.4: 20;
ntk2 = length(tk2);           % vp2, tk2:     upper source layer of upper crust 

tk3 = 5: 0.4: 15;             
ntk3 = length(tk3);           % vp2, tk3:     lower source layer of upper crust

vp4 = 6.1: 0.2: 8.3;
nvp4 = length(vp4);

tk4 = 10: 0.4: 25;
ntk4 = length(tk4);          % vp4, tk4:     lower crust

% vs1 = 1.0: 0.2: 3.4;
% nvs1 = length(vs1);
% 
% vs3 = 3.5: 0.2: 4.6;
% nvs3 = length(vs3);

%% 2. get the speed for P wave and thickness of each layer (7)
raypg = 0.1680: 0.0004: 0.1696;       % 1/5.87 ~= 0.1703,  max(raypg) < 1/max(vp2)
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
mul3pn = zeros(size3, 1);    %  pn
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

 
misfit = zeros(size1, size2, ntk3, size4);
tpn0pws = tp0 - raypn*498.6431;
% ii = 100;
% jj = 100;
% kk = 10;
% mm = 300;
for ii = 1: size1
    ii
    for jj = 1: size2
        
        % part 1, Pg(mod) ~= Pg(pick)
        distuse = [];
        pgmuse = [];
        pgduse = [];
        npg = zeros(nraypg, 1);
        temp = 0;
        for nn = 1: nraypg
            x = raypg(nn)*(div1pg(ii, nn)+div2pg(jj, nn));
            t = raypg(nn)*x+mul1pg(ii, nn)+mul2pg(jj, nn);
            [~, ind1]=min(abs(dist-(x-10)));
            [~, ind2]=min(abs(dist-(x+10)));
            pgmod = raypg(nn)*(dist(ind1: ind2)-x)+t;
            npg(nn) = ind2-ind1+1;
            temp = temp+sum((pgmod-pg(ind1: ind2)).^2)/npg(nn);
            distuse = [distuse; dist(ind1: ind2)];
            pgmuse = [pgmuse; pgmod];
            pgduse = [pgduse; pg(ind1: ind2)];            
        end
        misfit(ii, jj, :, :) = misfit(ii, jj, :, :) +temp/nraypg;
        
        % part 2, pPn-Pn(mod) ~= pPn-Pn(stack)
        dtpmod = 2*(mul1pn(ii)+mul2pn(jj));
        misfit(ii, jj, :, :) = misfit(ii, jj, :, :)+(dtpmod- dtp1).^2;
        
        [~, ivp2] = ind2sub([ntk2, nvp2], jj);
        for kk = 1: ntk3
            for mm = 1: size4                
                indice = (ivp2-1)*ntk3 + kk;
                
                % part 3, Pn intercept with x = 0
                xmin = raypn*(div1pn(ii)+div2pn(jj)+2*div3pn(indice)+2*div4pn(mm));
                tmin = raypn*xmin + mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                tpn0mod = mul1pn(ii)+mul2pn(jj)+2*mul3pn(indice)+2*mul4pn(mm);
                misfit(ii, jj, kk, mm) = misfit(ii, jj, kk, mm)+(tpn0mod-tpn0pws).^2;
                
                % part 4, Pn-Pg
                pnmuse = raypn*(distuse-xmin)+tmin;
                dtmod = pnmuse - pgmuse;
                pnduse = (distuse- 498.6431)/vpn + tp0;
                dtdata = pnduse -pgduse;
                dmd = dtmod - dtdata;
                temp = 0;
                for nn = 1:nraypg
                    inde = sum(npg(1:nn));
                    inds = inde-npg(nn)+1;
                    temp = temp+sum((dmd(inds: inde)).^2)/npg(nn);
                end
                misfit(ii, jj, kk, mm) = misfit(ii, jj, kk, mm)+temp/nraypg;
            end
        end
    end
end

[minval, ind] = min(misfit(:));
[is4, is3, is2, is1] = ind2sub([size4, ntk3, size2, size1], ind);
[itk1, ivp1] = ind2sub([ntk1, nvp1], is1);
[itk2, ivp2] = ind2sub([ntk2, nvp2], is2);
itk3 = is3;
[itk4, ivp4] = ind2sub([ntk4, nvp4], is4);
bestvp1 = vp1(ivp1);
bestvp2 = vp2(ivp2);
bestvp4 = vp4(ivp4);
besttk1 = tk1(itk1);
besttk2 = tk2(itk2);
besttk3 = tk3(itk3);
besttk4 = tk4(itk4);
bestdep = besttk1+besttk2;

save('p_wave_para.mat', 'besttk1', 'bestdep', 'bestvp1', 'besttk2', 'besttk3', 'bestvp4', 'bestvp2', 'besttk4', ...
         'ivp1', 'ivp2', 'ivp4', 'itk1', 'itk2', 'itk3', 'itk4', 'minval', 'ind', 'vp1', 'vp2', 'vp4', 'tk1','tk2', 'tk3', 'tk4', ...
         'is1', 'is2', 'is3', 'is4', 'size4', 'ntk3', 'size2', 'size1');


%% 3. get the speed for S wave of each layer (2), which is vs1 and vs3

% TWO ways to get vs1 and vs3

% %% 3.1     WAY 1: solve equations
% raysn =  1/vsn; 
% vs2 = 3.4;    
% % get vs1, use sPn-Pn
% etap1 = sqrt((1.0/bestvp1)^2-raypn^2);
% etap2 = sqrt((1.0/bestvp2)^2-raypn^2);
% %     etasp1 = sqrt((1.0/vs1)^2-rayp2^2);        % vs1 unknown
% etasp2 = sqrt((1.0/vs2)^2-raypn^2);
% temp1 = (besttk2) * (etap2+etasp2);
% temp2 = ( dtp2 - temp1 )/besttk1 - etap1; 
% vs1 = 1/sqrt( temp2^2 + raypn^2 );
% vs1 = roundn(vs1, -1);
%     
% % get vs3, use Sn intercept with x=0
% etas1 = sqrt((1.0/vs1)^2-raysn^2);
% etas2 = sqrt((1.0/vs2)^2-raysn^2);
% %     etas3 = sqrt((1.0/vs3)^2-rays2^2);     % vs3 unknown
% %     xsnmin = rays2 * ( besttk1/etas1 + (besttk1+2*besttk2-bestdep)/etas2 + 2*besttk3/etas3 );
% %     tsnmin = rays2 * xsnmin +besttk1*etas1 + (besttk1+2*besttk2-bestdep)*etas2 + 2*besttk3*etas3;
% %     tsn0mod = tsnmin - rays2*xsnmin;           % model determined pn curve's intercept point with time axes
% temp1 = ts0 - raysn*498.6431;             % PWS pn curve's intercept point with time axes
% temp2 = besttk1*etas1 + (besttk2+2*besttk3)*etas2;
% temp3 = (temp1 - temp2)/ (2*besttk4);
% vs3 = 1/sqrt( temp3^2 + raysn^2 );
% vs3 = roundn(vs3, -1);
%     
% save('s1_wave_para.mat', 'vs1', 'vs3');



%% 3.2      WAY 2:  min. misfit of model and data travel time curve
% set para. range 
vs1 = 1.5: 0.1: 3.0;
nvs1 = length(vs1);

vs2 = 3.25: 0.1: 3.45;
nvs2 = length(vs2);

vs4 = 3.6: 0.1: 4.6;
nvs4 = length(vs4);

raysn =  1/vsn;
raysg = 0.2844: 0.0005: 0.2864;    % 1/3.48 ~= 0.2874
nraysg = length(raysg);

% 3.2.1 sediment layer
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

% 3.2.2 upper crust
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

% 3.2.3 lower crust
% sn
div4sn = zeros(nvs4, 1);
mul4sn = zeros(nvs4, 1);      
for ii = 1:nvs4
    etasn = sqrt((1.0/vs4(ii))^2-raysn^2);
    div4sn(ii) = besttk4/etasn;
    mul4sn(ii) = besttk4*etasn; 
end


misfit = zeros(nvs1, nvs2, nvs4);
tsn0pws = ts0 - raysn*498.6431;             % PWS sn curve's intercept point with time axes
for ii = 1: nvs1
    ii
    for jj = 1: nvs2

        % part 1, Sg(mod) ~= Sg(pick)
        distuse = [];
        sgmuse = [];
        sgduse = [];
        nsg = zeros(nraysg, 1);
        temp = 0;
        for nn = 1: nraysg
            x = raysg(nn)*(div1sg(ii, nn)+div2sg(jj, nn));
            t = raysg(nn)*x+mul1sg(ii, nn)+mul2sg(jj, nn);
            [~, ind1]=min(abs(dist-(x-10)));
            [~, ind2]=min(abs(dist-(x+10)));
            sgmod = raysg(nn)*(dist(ind1: ind2)-x)+t;
            nsg(nn) = ind2-ind1+1;
            temp = temp + sum((sgmod-sg(ind1: ind2)).^2)/nsg(nn);
            distuse = [distuse; dist(ind1: ind2)];
            sgmuse = [sgmuse; sgmod];
            sgduse = [sgduse; sg(ind1: ind2)];            
        end
        misfit(ii, jj, :) = misfit(ii, jj, :) +temp/nraysg;
                
        for kk = 1:nvs4           
                
            % part 2, Sn intercept with x = 0
            xmin = raysn*(div1sn(ii)+div2sn(jj)+2*div3sn(jj)+2*div4sn(kk));
            tmin = raysn*xmin + mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            tsn0mod = mul1sn(ii)+mul2sn(jj)+2*mul3sn(jj)+2*mul4sn(kk);
            misfit(ii, jj, kk) = misfit(ii, jj, kk)+(tsn0mod-tsn0pws).^2;
                
            % part 3, Sn-Sg
            snmuse = raysn*(distuse-xmin)+tmin;
            dtmod = snmuse - sgmuse;
            snduse = (distuse - 498.6431)/vsn + ts0;
            dtdata = snduse -sgduse;
            dmd = dtmod - dtdata;
            temp = 0;
            for nn = 1:nraysg
                inde = sum(nsg(1:nn));
                inds = inde-nsg(nn)+1;
                temp = temp+sum((dmd(inds: inde)).^2)/nsg(nn);
            end
            misfit(ii, jj, kk) = misfit(ii, jj, kk)+temp/nraysg;
        end
    end
end

[minval, ind] = min(misfit(:));
[ivs4, ivs2, ivs1] = ind2sub([nvs4, nvs2, nvs1], ind);
bestvs1 = vs1(ivs1);
bestvs2 = vs2(ivs2);
bestvs4 = vs4(ivs4);

save('s2_wave_para.mat', 'bestvs1', 'bestvs2', 'bestvs4', 'ivs1', 'ivs2', 'ivs4', 'minval', 'ind', 'vs1', 'vs2', 'vs4'); 













