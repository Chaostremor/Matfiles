% travel time curve fitting
% 
% USAGE:
%     construct a four-layered model, sediment, upper crust, lower crust and mantle,
%     which the source is located in the upper crust
%     use the known and unknown parameters to calculate the theoratical curve and fit the real curve best
% Model:
%      tk1, vp1, vs1
%      tk2, vp2, vs2
%      tk3, vp3, vs3
%            vp4, vs4
% Author: C.Song,  2017.5.3

clear; 
%% 1. settle up the models
% load manually picked pg and sg arrival times named 'pg' and 'sg' 
load('pgsgpicktime.mat');
ndist = length(dist);

% results from PWS
vpn = 8.40;                          % vel. of Pn
tp0 = 72.12;                         % time of Pn at refdist 498.6431
dtp1 = 5.24;                         % difference time between Pn and pPn
dtp2 = 7.57;                         % difference time between Pn and sPn
vsn = 4.61;                           % vel. of Sn
ts0 = 122.8;                          % time of Sn at refdist 498.6431
dts1 = 8.0;                            % difference time between Sn and sSn

% known parameters
vp2 = 5.87;
vs2 = 3.48;
vp4 = vpn;
vs4 = vsn;

a1 = 5.969;
b1 = 35.05;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);           % pg travel time picks fitting curve
pncurve = (distcurve - 498.6431)/vp4 + tp0;    % PWS pn curve 

a2 = 10.34;
b2 = 35.95;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
sncurve = (distcurve - 498.6431)/vs4 + ts0;

figure
plot(pg, dist, 'k.', 'MarkerSize', 8); hold on;
plot(sg, dist, 'k.', 'MarkerSize', 8); hold on;
plot(pgcurve, distcurve, 'b-'); hold on;
plot(sgcurve, distcurve, 'b-'); hold on;
plot(pncurve, distcurve, 'r-'); hold on;
plot(sncurve, distcurve, 'r-'); hold on;

syms y x
S = solve(y==5.969*sqrt((x/35.05)^2+1), y==(x-498.6431)/8.4+72.12);
% S = solve('x-y=3', 'x+y=6', 'x', 'y');
xx = double(S.x);
tt = double(S.y);
crsx = roundn(xx(1), -1);
crst = roundn(tt(1), -1);

% unknown ones
vp1 = 1.5: 0.1: 5.8;
nvp1 = length(vp1);
tk1 = 0.1: 0.2: 5;
ntk1 = length(tk1);
%dep = 5: 0.2: 25;         % whatever the thickness of sediment is, source locates down to it


%% 2. get the depth, thickness and vp of sediment first 
misfit1 = zeros(ntk1, 101, nvp1);
misfit2 = zeros(ntk1, 101, nvp1);
nsample =10000;
rayp = linspace(0, 1.0/vp2-0.00005, nsample)';
icount = 0;
for ii = 1:ntk1                                 % thickness of sediment, tk
    dep = tk1(ii): 0.2 : tk1(ii)+20;       % force the source in the upper crust
    ndep = length(dep);
    for jj = 1:ndep                            % source depth, dep
        for kk = 1:nvp1                       % vp of sediment, vp1
            icount = icount +1;
            
            xp = zeros(nsample,1);
            tp = zeros(nsample,1);
            for mm = 1:nsample
                 etap1 = sqrt((1.0/vp1(kk))^2-rayp(mm)^2);
                 etap2 = sqrt((1.0/vp2)^2-rayp(mm)^2);
                 xp(mm) = rayp(mm) * ( tk1(ii)/etap1+(dep(jj)-tk1(ii))/etap2 );
                 tp(mm) = rayp(mm)*xp(mm)+tk1(ii)*etap1+(dep(jj)-tk1(ii))*etap2;
            end
            pgfit = interp1(xp, tp, distcurve, 'spline');       % model determined pg curve 
            misfit1(ii, jj, kk) = sum((pgfit-pgcurve).^2);          
            
            dtpmod = 2*(dep(jj)-tk1(ii))*sqrt((1/vp2)^2-(1/vp4)^2) + 2*tk1(ii)*sqrt((1/vp2)^2-(1/vp4)^2);  % model determined difference between Pn and pPn
            misfit2(ii, jj, kk) = ( dtpmod- dtp1)^2;       % PWS difference
            
        end
    end
end
misfit = misfit1+ misfit2;
[~, index] = min(misfit(:));
[itk1, idep, ivp1] = ind2sub([ntk1, ndep, nvp1], index);
besttk1 = tk1(itk1);
bestdep =dep(idep);
bestvp1 = vp1(ivp1); 

% now remain tk2, tk3, vp3

%% 3. get thickness of upper crust tk2, thickness and vp of lower crust tk3, vp3
% may be a underdetermined problem 
% unknown ones
vp3 = 5.9: 0.1: 8.3;
nvp3 = length(vp3);
tk2 = bestdep-besttk1: 0.2: 30;
ntk2 = length(tk2);
tk3 = 10: 0.2: 30;
ntk3 = length(tk3);
misfit1 = zeros(ntk2, ntk3, nvp3);
misfit2 = zeros(ntk2, ntk3, nvp3);
rayp = 1/vp4;
icount = 0;
for ii = 1:ntk2
    for jj = 1:ntk3
        for kk = 1:nvp3
            icount = icount +1;
            
            etap1 = sqrt((1.0/bestvp1)^2-rayp^2);
            etap2 = sqrt((1.0/vp2)^2-rayp^2);
            etap3 = sqrt((1.0/vp3(kk))^2-rayp^2);
            xpnmin = rayp * ( besttk1/etap1 + (besttk1+2*tk2(ii)-bestdep)/etap2 + 2*tk3(jj)/etap3 );
            tpnmin = rayp * xpnmin +besttk1*etap1 + (besttk1+2*tk2(ii)-bestdep)*etap2 + 2*tk3(jj)*etap3;
            tpn0mod = tpnmin - rayp*xpnmin;           % model determined pn curve's intercept point with time axes
            tpn0pws = 72.12 - rayp*498.6431;             % PWS pn curve's intercept point with time axes
            misfit1(ii, jj, kk) = (tpn0mod - tpn0pws)^2;
            
            tpn = rayp * crsx +besttk1*etap1 + (besttk1+2*tk2(ii)-bestdep)*etap2 + 2*tk3(jj)*etap3;  % model determined pn's intercept with pg
            misfit2(ii, jj, kk) = (tpn - crst)^2;          % PWS pn curve's intercept with pg
        end
    end
end

misfit = misfit1+ misfit2;
[~, index] = min(misfit(:));
[itk2, itk3, ivp3] = ind2sub([ntk2, ntk3, nvp3], index);
besttk2 = tk2(itk2);
besttk3 =tk3(itk3);
bestvp3 = vp3(ivp3); 

save('travel1.mat');









