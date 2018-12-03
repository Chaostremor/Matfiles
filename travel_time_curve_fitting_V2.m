% travel time curve fitting V2
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
%
%
% Author: C.Song,  2017.5.7

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
vp2 = 5.87;
vs2 = 3.48;
vp4 = vpn;
vs4 = vsn;

a1 = 5.969;
b1 = 35.05;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);           % pg travel time picks fitting curve
pncurve = (distcurve - 498.6431)/vp4 + tp0;    % PWS pn curve
pn = (dist- 498.6431)/vp4 + tp0;

a2 = 10.34;
b2 = 35.95;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
sncurve = (distcurve - 498.6431)/vs4 + ts0;
sn = (dist - 498.6431)/vs4 + ts0;

% figure
% plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;
% plot(sg, dist, 'k.', 'MarkerSize', 8); hold on;
% plot(pgcurve, distcurve, 'b-'); hold on;
% plot(sgcurve, distcurve, 'b-'); hold on;
% plot(pncurve, distcurve, 'r-'); hold on;
% plot(sncurve, distcurve, 'r-'); hold on;

syms y x
S = solve(y==5.969*sqrt((x/35.05)^2+1), y==(x-498.6431)/8.4+72.12);
% S = solve('x-y=3', 'x+y=6', 'x', 'y');
xx = double(S.x);
tt = double(S.y);
crsx = roundn(xx(1), -1);
crst = roundn(tt(1), -1);

% unknown ones  (8)
vp1 = 2.5: 0.2: 5;
nvp1 = length(vp1);

vp3 = 6: 0.2: 8.3;
nvp3 = length(vp3);

tk1 = 0.1: 0.2: 5;
ntk1 = length(tk1);

tk2 = 12.5: 0.4: 30;
ntk2 = length(tk2);

tk3 = 12.5: 0.4: 30;      % tk2 and tk3 make sure that source locates up to Moho
ntk3 = length(tk3);

dep = 5: 0.2: 25;         % whatever the thickness of sediment is, source locates down to it
ndep = length(dep);

% vs1 = 1.0: 0.2: 3.4;
% nvs1 = length(vs1);
% 
% vs3 = 3.5: 0.2: 4.6;
% nvs3 = length(vs3);


%% 2. get the speed for P wave and thickness of each layer (6) 

nsample2 = 1000;       % for Pn
rayp2 =  1/vp4;
icount = 1;
size = ntk1*ndep*nvp1*ntk2*ntk3*nvp3;
misfit = zeros(size, 1);
% ii =10;
% jj = 10;
% kk =10;
% ll = 10;
% mm = 10;
% nn = 10;
for ii = 1:ntk1                                 % thickness of sediment, tk1
    for jj = 1:ndep                            % source depth, dep
        for kk = 1:nvp1                      % vp of sediment, vp1            
            for ll = 1:ntk2                     % thickness of upper crust, tk2
                for mm = 1:ntk3             % thickness of lower crust, tk3
                    for nn = 1:nvp3          % vp of lower crust, vp3
                        if ( dep(jj) <= tk1(ii) + tk2(ll) )         % source is in the upper crust
                            iflag =1;
                            
                            % part 1, Pg(mod) ~= Pg(fit)
                            nsample1 = 10000;     % for Pg 
                            rayp1 = linspace(0, 1.0/vp2-0.00005, nsample1)';
                            xp = zeros(nsample1,1);
                            tp = zeros(nsample1,1);
                            for qq = 1:nsample1
                                etap1 = sqrt((1.0/vp1(kk))^2-rayp1(qq)^2);
                                etap2 = sqrt((1.0/vp2)^2-rayp1(qq)^2);
                                xp(qq) = rayp1(qq) * ( tk1(ii)/etap1 + (dep(jj)-tk1(ii))/etap2 );
                                tp(qq) = rayp1(qq)*xp(qq) + tk1(ii)*etap1 + (dep(jj)-tk1(ii))*etap2;                              
                            end
                            pgfit = interp1(xp, tp, dist, 'spline');       % model determined pg curve
                            temp1 = sum((pgfit-pg).^2)/ndist; 
                            
                            % part 2, pPn-Pn(mod) ~= pPn-Pn(stack)
                            etap1 = sqrt((1.0/vp1(kk))^2-rayp2^2);
                            etap2 = sqrt((1.0/vp2)^2-rayp2^2);
                            dtpmod = 2*tk1(ii)*etap1 + 2*(dep(jj)-tk1(ii))*etap2;  % model determined difference between Pn and pPn
                            temp2 = ( dtpmod- dtp1)^2;       % PWS difference                            
                            
                            % part 3, Pn intercept with x = 0                            
%                             etap1 = sqrt((1.0/vp1(kk))^2-rayp2^2);
%                             etap2 = sqrt((1.0/vp2)^2-rayp2^2);
                            etap3 = sqrt((1.0/vp3(nn))^2-rayp2^2);
                            xpnmin = rayp2 * ( tk1(ii)/etap1 + (tk1(ii)+2*tk2(ll)-dep(jj))/etap2 + 2*tk3(mm)/etap3 );
                            tpnmin = rayp2 * xpnmin +tk1(ii)*etap1 + (tk1(ii)+2*tk2(ll)-dep(jj))*etap2 + 2*tk3(mm)*etap3;
                            tpn0mod = tpnmin - rayp2*xpnmin;           % model determined pn curve's intercept point with time axes
                            tpn0pws = tp0 - rayp2*498.6431;             % PWS pn curve's intercept point with time axes
                            temp3 = (tpn0mod - tpn0pws)^2;
                            
                            % part 4, Pn-Pg
                            dxpn = linspace(0, 400, nsample2)';
                            xpn = xpnmin +dxpn;
                            tpn = tpnmin +rayp2 .*dxpn;
                            pnfit = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
                            dtmodel = pnfit - pgfit;
                            dtdata = pn - pg;
                            temp4 = sum((dtmodel - dtdata).^2)/ndist;
                            
                            % sum all constraint parts
                            misfit(icount) = temp1+ temp2 +temp3+temp4;
                            icount = icount +1;
                            %fprintf('%d', icount);
                            
                        elseif  ( dep(jj) > tk1(ii) + tk2(ll) && dep(jj) < tk1(ii) + tk2(ll) +tk3(mm) )   % source is in the lower crust                                                         xp = zeros(nsample1,1);
                            iflag =2;
                            
                            % part 1, Pg(mod) ~= Pg(fit)
                            nsample1 = 10000;     % for Pg
                            rayp1 = linspace(0, 1.0/vp3(nn)-0.00005, nsample1)';
                            xp = zeros(nsample1,1);
                            tp = zeros(nsample1,1);
                            for qq = 1:nsample1
                                etap1 = sqrt((1.0/vp1(kk))^2-rayp1(qq)^2);
                                etap2 = sqrt((1.0/vp2)^2-rayp1(qq)^2);
                                etap3 = sqrt((1.0/vp3(nn))^2-rayp1(qq)^2);
                                xp(qq) = rayp1(qq) * ( tk1(ii)/etap1 + tk2(ll)/etap2 + (dep(jj)-tk1(ii)-tk2(ll))/etap3 );
                                tp(qq) = rayp1(qq)*xp(qq) + tk1(ii)*etap1 + tk2(ll)*etap2 + (dep(jj)-tk1(ii)-tk2(ll))*etap3;                             
                            end
                            pgfit = interp1(xp, tp, dist, 'spline');       % model determined pg curve
                            temp1 = sum((pgfit-pg).^2)/ndist;
                            
                            % part 2, pPn-Pn(mod) ~= pPn-Pn(stack)
                            etap1 = sqrt((1.0/vp1(kk))^2-rayp2^2);
                            etap2 = sqrt((1.0/vp2)^2-rayp2^2);
                            etap3 = sqrt((1.0/vp3(nn))^2-rayp2^2);
                            dtpmod = 2*tk1(ii)*etap1 + 2*tk2(ll)*etap2 + 2*(dep(jj)-tk1(ii)-tk2(ll))*etap3;  % model determined difference between Pn and pPn
                            temp2 = ( dtpmod- dtp1)^2;       % PWS difference                            
                            
                            % part 3, Pn intercept with x = 0                            
                            xpnmin = rayp2 * ( tk1(ii)/etap1 + tk2(ll)/etap2 + (tk1(ii)+tk2(ll)+2*tk3(mm)-dep(jj))/etap3 );
                            tpnmin = rayp2 * xpnmin + tk1(ii)*etap1 + tk2(ll)*etap2 + (tk1(ii)+tk2(ll)+2*tk3(mm)-dep(jj))*etap3;
                            tpn0mod = tpnmin - rayp2*xpnmin;           % model determined pn curve's intercept point with time axes
                            tpn0pws = tp0 - rayp2*498.6431;             % PWS pn curve's intercept point with time axes
                            temp3 = (tpn0mod - tpn0pws)^2;
                            
                            % part 4, Pn-Pg
                            dxpn = linspace(0, 400, nsample2)';
                            xpn = xpnmin +dxpn;
                            tpn = tpnmin +rayp2 .*dxpn;
                            pnfit = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
                            dtmodel = pnfit - pgfit;
                            dtdata = pn - pg;
                            temp4 = sum((dtmodel - dtdata).^2)/ndist;
                            
                            % sum all constraint parts
                            misfit(icount) = temp1+ temp2 +temp3+temp4;
                            icount = icount +1;
                            %fprintf('%d', icount);
                            
                        end
                    end
                end
            end
        end
    end
end

[minval, index] = min(misfit(:));
[ivp3, itk3, itk2, ivp1, idep, itk1] = ind2sub([nvp3, ntk3, ntk2, nvp1, ndep, ntk1], index);
besttk1 = tk1(itk1);
bestdep =dep(idep);
bestvp1 = vp1(ivp1);
besttk2 = tk2(itk2);
besttk3 = tk3(itk3);
bestvp3 = vp3(ivp3);

save('p_wave_para.mat', 'besttk1', 'bestdep', 'bestvp1', 'besttk2', 'besttk3', 'bestvp3', 'ivp3', 'itk3', 'itk2', 'ivp1', 'idep', 'itk1', ...
        'minval', 'index', 'tk1', 'dep', 'vp1', 'tk2', 'tk3', 'vp3');

%% 3. get the speed for S wave of each layer (2), which is vs1 and vs3

% TWO ways to get vs1 and vs3

% %% 3.1     WAY 1: solve equations
% rays2 =  1/vs4; 
% if ( bestdep <= besttk1 + besttk2 )    % if source is in the upper crust
%     
%     % get vs1, use sPn-Pn
%     etap1 = sqrt((1.0/bestvp1)^2-rayp2^2);
%     etap2 = sqrt((1.0/vp2)^2-rayp2^2);
% %     etasp1 = sqrt((1.0/vs1)^2-rayp2^2);        % vs1 unknown
%     etasp2 = sqrt((1.0/vs2)^2-rayp2^2);
%     temp1 = (bestdep-besttk1) * (etap2+etasp2);
%     temp2 = ( dtp2 - temp1 )/besttk1 - etap1; 
%     vs1 = 1/sqrt( temp2^2 + rayp2^2 );
%     vs1 = roundn(vs1, -1);
%     
%     % get vs3, use Sn intercept with x=0
%     etas1 = sqrt((1.0/vs1)^2-rays2^2);
%     etas2 = sqrt((1.0/vs2)^2-rays2^2);
% %     etas3 = sqrt((1.0/vs3)^2-rays2^2);     % vs3 unknown
% %     xsnmin = rays2 * ( besttk1/etas1 + (besttk1+2*besttk2-bestdep)/etas2 + 2*besttk3/etas3 );
% %     tsnmin = rays2 * xsnmin +besttk1*etas1 + (besttk1+2*besttk2-bestdep)*etas2 + 2*besttk3*etas3;
% %     tsn0mod = tsnmin - rays2*xsnmin;           % model determined pn curve's intercept point with time axes
%     temp1 = ts0 - rays2*498.6431;             % PWS pn curve's intercept point with time axes
%     temp2 = besttk1*etas1 + (besttk1+2*besttk2-bestdep)*etas2;
%     temp3 = (temp1 - temp2)/ (2*besttk3);
%     vs3 = 1/sqrt( temp3^2 + rays2^2 );
%     vs3 = roundn(vs3, -1);
%     
% elseif ( bestdep > besttk1 + besttk2 && bestdep < besttk1 + besttk2 + besttk3 )     % source is in the lower crust 
%     
%     % use sPn-Pn and Sn intercept with x=0 together
%     syms vs1 vs3
%     etap1 = sqrt((1.0/bestvp1)^2-rayp2^2);
%     etap2 = sqrt((1.0/vp2)^2-rayp2^2);
%     etap3 = sqrt((1.0/bestvp3)^2-rayp2^2);
%     etasp1 = sqrt((1.0/vs1)^2-rayp2^2);        % vs1 unknown
%     etasp2 = sqrt((1.0/vs2)^2-rayp2^2);
%     etasp3 = sqrt((1.0/vs3)^2-rayp2^2);        % vs3 unknown   
%     S1 = besttk1*(etap1+etasp1) + besttk2*(etap2+etasp2) + (bestdep-besttk1-besttk2)*(etap3+etasp3) - dtp2;
%     etas1 = sqrt((1.0/vs1)^2-rays2^2);      % vs1 unknown
%     etas2 = sqrt((1.0/vs2)^2-rays2^2);
%     etas3 = sqrt((1.0/vs3)^2-rays2^2);     % vs3 unknown
%     S2 = besttk1*etas1 + besttk2*etas2 + (besttk1+besttk2+2*besttk3-bestdep)*etas3 - ts0;
%     S = solve(S1, S2);
%     vs1 = roundn( double(S.vs1), -1 );
%     vs3 = roundn( double(S.vs3), -1);
%     
% end
% save('s_wave_para.mat', 'vs1', 'vs3');

%% 3.2      WAY 2:  min. misfit of model and data travel time curve
% set para. range 
vs1 = 1.0: 0.2: 3.4;
nvs1 = length(vs1);
    
vs3 = 3.5: 0.2: 4.6;
nvs3 = length(vs3);

nsample2 = 1000;       % for Sn
rays2 =  1/vs4;
icount = 1;
size = nvs1*nvs3;
misfit = zeros(size, 1);
if ( bestdep <= besttk1 + besttk2 )    % if source is in the upper crust
    for ii = 1:nvs1
        for jj = 1:nvs3
            
            % part 1, Sg(mod) ~= Sg(fit)
            nsample1 = 10000;     % for Sg
            rays1 = linspace(0, 1.0/vs2-0.0000005, nsample1)';
            xs = zeros(nsample1,1);
            ts = zeros(nsample1,1);
            for qq = 1:nsample1
                etas1 = sqrt((1.0/vs1(ii))^2-rays1(qq)^2);
                etas2 = sqrt((1.0/vs2)^2-rays1(qq)^2);
                xs(qq) = rays1(qq) * ( besttk1/etas1 + (bestdep-besttk1)/etas2 );
                ts(qq) = rays1(qq)*xs(qq) + besttk1*etas1 + (bestdep-besttk1)*etas2;                              
            end
            sgfit = interp1(xs, ts, dist, 'spline');       % model determined sg curve
            temp1 = sum((sgfit-sg).^2)/ndist; 
            
            % part 2 , Sn intercept with x = 0                            
            etas1 = sqrt((1.0/vs1(ii))^2-rays2^2);
            etas2 = sqrt((1.0/vs2)^2-rays2^2);
            etas3 = sqrt((1.0/vs3(jj))^2-rays2^2);
            xsnmin = rays2 * ( besttk1/etas1 + (besttk1+2*besttk2-bestdep)/etas2 + 2*besttk3/etas3 );
            tsnmin = rays2 * xsnmin +besttk1*etas1 + (besttk1+2*besttk2-bestdep)*etas2 + 2*besttk3*etas3;
            tsn0mod = tsnmin - rays2*xsnmin;           % model determined sn curve's intercept point with time axes
            tsn0pws = ts0 - rays2*498.6431;             % PWS sn curve's intercept point with time axes
            temp2 = (tsn0mod - tsn0pws)^2;
            
            % part 3, Sn-Sg
            dxsn = linspace(0, 400, nsample2)';
            xsn = xsnmin +dxsn;
            tsn = tsnmin +rays2 .*dxsn;           
            snfit = interp1(xsn, tsn, dist, 'spline');       % model determined sn curve
            dtmodel = snfit - sgfit;
            dtdata = sn - sg;
            temp3 = sum((dtmodel - dtdata).^2)/ndist;
            
            % sum all constraint parts
            misfit(icount) = temp1+ temp2 + temp3;
            icount = icount +1;
            %fprintf('%d', icount);
        end
    end
    
elseif ( bestdep > besttk1 + besttk2 && bestdep < besttk1 + besttk2 + besttk3 )     % source is in the lower crust
    for ii = 1:nvs1
        for jj = 1:nvs3
             
            % part 1, Sg(mod) ~= Sg(fit)
            nsample1 = 10000;     % for Sg
            rays1 = linspace(0, 1.0/vs3(jj)-0.0000005, nsample1)';
            xs = zeros(nsample1,1);
            ts = zeros(nsample1,1);
            for qq = 1:nsample1
                etas1 = sqrt((1.0/vs1(ii))^2-rays1(qq)^2);
                etas2 = sqrt((1.0/vs2)^2-rays1(qq)^2);
                etas3 = sqrt((1.0/vs3(jj))^2-rays1(qq)^2);
                xs(qq) = rays1(qq) * ( besttk1/etas1 + besttk2/etas2 + (bestdep-besttk1-besttk2)/etas3 );
                ts(qq) = rays1(qq)*xs(qq) + besttk1*etas1 + besttk2*etas2 + (bestdep-besttk1-besttk2)*etas3;                             
            end
            sgfit = interp1(xs, ts, dist, 'spline');       % model determined pg curve
            temp1 = sum((sgfit-sg).^2)/ndist;
                                                                        
            % part 2, Sn intercept with x = 0
            etas1 = sqrt((1.0/vs1(ii))^2-rays2^2);
            etas2 = sqrt((1.0/vs2)^2-rays2^2);
            etas3 = sqrt((1.0/vs3(jj))^2-rays2^2);
            xsnmin = rays2 * ( besttk1/etas1 + besttk2/etas2 + (besttk1+besttk2+2*besttk3-bestdep)/etas3 );
            tsnmin = rays2 * xsnmin + besttk1*etas1 + besttk2*etas2 + (besttk1+besttk2+2*besttk3-bestdep)*etas3;
            tsn0mod = tsnmin - rays2*xsnmin;           % model determined Sn curve's intercept point with time axes
            tsn0pws = ts0 - rays2*498.6431;             % PWS Sn curve's intercept point with time axes
            temp2 = (tsn0mod - tsn0pws)^2;
                            
            % part 3, Sn-Sg
            dxsn = linspace(0, 400, nsample2)';
            xsn = xsnmin +dxsn;
            tsn = tsnmin +rays2 .*dxsn;
            snfit = interp1(xsn, tsn, dist, 'spline');       % model determined Sn curve
            dtmodel = snfit - sgfit;
            dtdata = sn - sg;
            temp3 = sum((dtmodel - dtdata).^2)/ndist;
                            
            % sum all constraint parts
            misfit(icount) = temp1+ temp2 + temp3;
            icount = icount +1;
            %fprintf('%d', icount);
        end
    end
    
end

[minval, index] = min(misfit(:));
[ivs3, ivs1] = ind2sub([nvs3, nvs1], index);
bestvs1 = vs1(ivs1);
bestvs3 = vs3(ivs3);

save('s_wave_para.mat', 'bestvs1', 'bestvs3', 'ivs1', 'ivs3', 'minval', 'index', 'vs1', 'vs3'); 
