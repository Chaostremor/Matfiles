% velocity structure inversion v1.0
%
% USAGE:
%     use least-square method to invert velocity structure based on an initial
%     model
% Author: C. Song, 2017.5.10
%


%% 1. read data
datadir = 'G:\Alxa\nodecimate\3test600\' ;    % 数据所在目录
fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
[sa, sb] = size(stnm);
% test to get size of data
ii=1;
filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % 获取文件路径
[ttest, datatest, SAChdrtest] = fget_sac(filename) ;
disttest = SAChdrtest.evsta.dist ;   % 注意头段变量的存储方式
dt = SAChdrtest.times.delta ;
npts = SAChdrtest.data.trcLen ; 
%
[sc, sd] = size(ttest);
t = zeros(sc,sa);          % t是时间，data是数据
data = zeros(sc,sa);
dist = zeros(sa,1);       % dist是头段中的震中距
for ii = 1:sa
    filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % 获取文件路径
    [t(:,ii), data(:,ii), SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    dist(ii) = SAChdr.evsta.dist ;   % 注意头段变量的存储方式
end

% initial model
vs = [3.04; 3.51; 3.77; 4.756];
vp = [5.09208; 5.90612; 6.35644; 8.06419];
tk= [5.1; 20.0; 15.0];
dep = 16.0;

%% 2. Pg

% for example, as for Pg wave
% assume source is in upper crust of three-layed model
% Pg wave traveltime curve
nsample1 = 10000;      
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for qq = 1:nsample1    
    etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
    xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
    tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
pg0 = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
% rayp = zeros(sa, 1);
% syms x
% etap1 = sqrt((1.0/vp(1))^2-x^2);
% etap2 = sqrt((1.0/vp(2))^2-x^2);
% for i = 1:sa
% S1 = dist(i) - x * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
% S = solve(S1);
% rayp(i) = double(S.x);
% end

% tpg = f(h, tk(1), vp(1)), use total differential equation to get dtpg
ddep = 0.1;
dtk1 = 0.1;
dvp1 = 0.1;
dvp2 = 0.1;
% drayp = 0.001;

% for each distance
% compute partial derivative
% part 1: dtdh
dep = dep+ddep;   
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for qq = 1:nsample1    
    etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
    xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
    tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
pgdep = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
dtpgdh = (pgdep - pg0)./ddep;
dep = dep-ddep;

% part 2: dtdtk1
tk(1) = tk(1)+dtk1;    
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for qq = 1:nsample1    
    etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
    xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
    tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
pgtk1 = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
dtpgdtk1 = (pgtk1 - pg0)./dtk1;
tk(1) = tk(1)-dtk1;

% part 3: dtdvp1
vp(1) = vp(1)+dvp1;      
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for qq = 1:nsample1    
    etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
    xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
    tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
pgvp1 = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
dtpgdvp1 = (pgvp1 - pg0)./dvp1;
vp(1) = vp(1)-dvp1;

% part 4: dtdvp2
vp(2) = vp(2)+dvp2;      
rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for qq = 1:nsample1    
    etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
    xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
    tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
pgvp2 = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
dtpgdvp2 = (pgvp2 - pg0)./dvp2;
vp(2) = vp(2)-dvp2;

% % part 5: dtdrayp
% vp(2) = vp(2)-dvp2;
% nsample1 = 10000;      
% rayp1 = linspace(0, 1.0/vp(2)-0.00005, nsample1)';
% xpg = zeros(nsample1,1);
% tpg = zeros(nsample1,1);
% for qq = 1:nsample1    
%     etap1 = sqrt((1.0/vp(1))^2-rayp1(qq)^2);
%     etap2 = sqrt((1.0/vp(2))^2-rayp1(qq)^2);
%     xpg(qq) = rayp1(qq) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );  
%     tpg(qq) = rayp1(qq)*xpg(qq) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
% end
% pgp = interp1(xpg, tpg, dist, 'spline');       % model determined pg curve
% dtdtk1 = (pgtk1 - pg0)./dtk1;

% take all parts into account
% M = [deldep; deltk1; delvp1; delvp2];
G1 = [dtpgdh dtpgdtk1 zeros(sa, 2) dtpgdvp1 dtpgdvp2 zeros(sa, 6)];
D1 = dtpg;  
% dtpg comes from cross-correlation of phase Pg of data and synthetics

% LS solution: GM=D , M=((G'G)^-1)G'D


%% 3. Pn

rayp2 = 1/vp(4);
nsample2 = 1000;
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
dxpn = linspace(0, 400, nsample2)';
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pn0 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve

% tpn = f(h, tk(1~3), vp(1~4)), use total differential equation to get dtpn
ddep = 0.1;
dtk1 = 0.1;
dtk2 = 0.1;
dtk3 = 0.1;
dvp1 = 0.1;
dvp2 = 0.1;
dvp3 = 0.1;
dvp4 = 0.1;

% for each distance
% compute partial derivative
% part 1: dtpndh
dep = dep+ddep;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndep = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndh = (pndep-pn0)./ddep;
dep = dep-ddep;

% part 2: dtpndtk1
tk(1) = tk(1)+dtk1;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndtk1 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndtk1 = (pndtk1-pn0)./dtk1;
tk(1) = tk(1)-dtk1;

% part 3: dtpndtk2
tk(2) = tk(2)+dtk2;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndtk2 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndtk2 = (pndtk2-pn0)./dtk2;
tk(2) = tk(2)-dtk2;

% part 4: dtpndtk3
tk(3) = tk(3)+dtk3;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndtk3 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndtk3 = (pndtk3-pn0)./dtk3;
tk(3) = tk(3)-dtk3;

% part 5: dtpndvp1
vp(1) =vp(1)+dvp1;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndvp1 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndvp1 = (pndvp1-pn0)./dvp1;
vp(1) =vp(1)-dvp1;

% part 6: dtpndvp2
vp(2) =vp(2)+dvp2;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndvp2 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndvp2 = (pndvp2-pn0)./dvp2;
vp(2) =vp(2)-dvp2;

% part 7: dtpndvp3
vp(3) =vp(3)+dvp3;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndvp3 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndvp3 = (pndvp3-pn0)./dvp3;
vp(3) =vp(3)-dvp3;

% part 8: dtpndvp4
vp(4) =vp(4)+dvp4;
rayp2 = 1/vp(4);
etap1 = sqrt((1.0/vp(1))^2-rayp2^2);
etap2 = sqrt((1.0/vp(2))^2-rayp2^2);
etap3 = sqrt((1.0/vp(3))^2-rayp2^2);
xpnmin = rayp2 * ( tk(1)/etap1 + (tk(1)+2*tk(2)-dep)/etap2 + 2*tk(3)/etap3 );
tpnmin = rayp2 * xpnmin +tk(1)*etap1 + (tk(1)+2*tk(2)-dep)*etap2 + 2*tk(3)*etap3;
xpn = xpnmin +dxpn;
tpn = tpnmin +rayp2 .*dxpn;
pndvp4 = interp1(xpn, tpn, dist, 'spline');       % model determined pn curve
dtpndvp4 = (pndvp4-pn0)./dvp4;
vp(4) =vp(4)-dvp4;

% take all parts into account
G2 = [dtpndh dtpndtk1 dtpndtk2 dtpndtk3 dtpndvp1 dtpndvp2 dtpndvp3 dtpndvp4 zeros(sa, 4)];
D2 = dtpn;
% dtpn also comes from the cross-correlation of phase Pn of data and synthetics


%% 4. Sg     
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for qq = 1:nsample1    
    etas1 = sqrt((1.0/vs(1))^2-rays1(qq)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(qq)^2);
    xsg(qq) = rays1(qq) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );  
    tsg(qq) = rays1(qq)*xsg(qq) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
sg0 = interp1(xsg, tsg, dist, 'spline');       % model determined sg curve


% tsg = f(h, tk(1), vs(1)), use total differential equation to get dtsg
ddep = 0.1;
dtk1 = 0.1;
dvs1 = 0.1;
dvs2 = 0.1;

% for each distance
% compute partial derivative
% part 1: dtdh
dep = dep+ddep;    
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for qq = 1:nsample1    
    etas1 = sqrt((1.0/vs(1))^2-rays1(qq)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(qq)^2);
    xsg(qq) = rays1(qq) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );  
    tsg(qq) = rays1(qq)*xsg(qq) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
sgdep = interp1(xsg, tsg, dist, 'spline');       % model determined pg curve
dtsgdh = (sgdep - sg0)./ddep;
dep = dep-ddep;

% part 2: dtdtk1
tk(1) = tk(1)+dtk1;    
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for qq = 1:nsample1    
    etas1 = sqrt((1.0/vs(1))^2-rays1(qq)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(qq)^2);
    xsg(qq) = rays1(qq) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );  
    tsg(qq) = rays1(qq)*xsg(qq) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
sgtk1 = interp1(xsg, tsg, dist, 'spline');       % model determined pg curve
dtsgdtk1 = (sgtk1 - sg0)./dtk1;
tk(1) = tk(1)-dtk1;

% part 3: dtdvs1
vs(1) = vs(1)+dvs1;     
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for qq = 1:nsample1    
    etas1 = sqrt((1.0/vs(1))^2-rays1(qq)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(qq)^2);
    xsg(qq) = rays1(qq) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );  
    tsg(qq) = rays1(qq)*xsg(qq) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
sgvs1 = interp1(xsg, tsg, dist, 'spline');       % model determined pg curve
dtsgdvs1 = (sgvs1 - sg0)./dvs1;
vs(1) = vs(1)-dvs1;

% part 4: dtdvs2
vs(2) = vs(2)+dvs2;   
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for qq = 1:nsample1    
    etas1 = sqrt((1.0/vs(1))^2-rays1(qq)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(qq)^2);
    xsg(qq) = rays1(qq) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );  
    tsg(qq) = rays1(qq)*xsg(qq) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
sgvs2 = interp1(xsg, tsg, dist, 'spline');       % model determined pg curve
dtsgdvs2 = (sgvs2 - sg0)./dvs2;
vs(2) = vs(2)-dvs2;


% take all parts into account
% M = [deldep; deltk1; delvp1; delvp2];
G3 = [dtsgdh dtsgdtk1 zeros(sa, 6) dtsgdvs1 dtsgdvs2 zeros(sa, 2)];
D3 = dtsg;  
% dtsg comes from cross-correlation of phase Sg of data and synthetics

%% 6. Sn

rays2 = 1/vs(4);
nsample2 = 1000;
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
dxsn = linspace(0, 400, nsample2)';
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sn0 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve

% tpn = f(h, tk(1~3), vs(1~4)), use total differential equation to get dtpn
ddep = 0.1;
dtk1 = 0.1;
dtk2 = 0.1;
dtk3 = 0.1;
dvs1 = 0.1;
dvs2 = 0.1;
dvs3 = 0.1;
dvs4 = 0.1;

% for each distance
% compute partial derivative
% part 1: dtsndh
dep = dep+ddep;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndep = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndh = (sndep-sn0)./ddep;
dep = dep-ddep;

% part 2: dtsndtk1
tk(1) = tk(1)+dtk1;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndtk1 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndtk1 = (sndtk1-sn0)./dtk1;
tk(1) = tk(1)-dtk1;

% part 3: dtsndtk2
tk(2) = tk(2)+dtk2;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndtk2 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndtk2 = (sndtk2-sn0)./dtk2;
tk(2) = tk(2)-dtk2;

% part 4: dtsndtk3
tk(3) = tk(3)+dtk3;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndtk3 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndtk3 = (sndtk3-sn0)./dtk3;
tk(3) = tk(3)-dtk3;

% part 5: dtsndvs1
vs(1) =vs(1)+dvs1;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndvs1 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndvs1 = (sndvs1-sn0)./dvs1;
vs(1) =vs(1)-dvs1;

% part 6: dtsndvs2
vs(2) =vs(2)+dvs2;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndvs2 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndvs2 = (sndvs2-sn0)./dvs2;
vs(2) =vs(2)-dvs2;

% part 7: dtsndvs3
vs(3) =vs(3)+dvs3;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndvs3 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndvs3 = (sndvs3-sn0)./dvs3;
vs(3) =vs(3)-dvs3;

% part 8: dtsndvs4
vs(4) =vs(4)+dvs4;
rays2 = 1/vs(4);
etas1 = sqrt((1.0/vs(1))^2-rays2^2);
etas2 = sqrt((1.0/vs(2))^2-rays2^2);
etas3 = sqrt((1.0/vs(3))^2-rays2^2);
xsnmin = rays2 * ( tk(1)/etas1 + (tk(1)+2*tk(2)-dep)/etas2 + 2*tk(3)/etas3 );
tsnmin = rays2 * xsnmin +tk(1)*etas1 + (tk(1)+2*tk(2)-dep)*etas2 + 2*tk(3)*etas3;
xsn = xsnmin +dxsn;
tsn = tsnmin +rays2 .*dxsn;
sndvs4 = interp1(xsn, tsn, dist, 'spline');       % model determined pn curve
dtsndvs4 = (sndvs4-sn0)./dvs4;
vs(4) =vs(4)-dvs4;

% take all parts into account
G4 = [dtpndh dtpndtk1 dtpndtk2 dtpndtk3  zeros(sa, 4) dtsndvs1 dtsndvs2 dtsndvs3 dtsndvs4];
D2 = dtsn;
% dtsn also comes from the cross-correlation of phase Sn of data and synthetics


D = [D1; D2; D3; D4];
G = [G1; G2; G3; G4];
M =  (G'*G)\G'*D;
dep = dep+M(1);
tk(1) = tk(1)+M(2);
tk(2) = tk(2)+M(3);
tk(3) = tk(3)+M(4);
vp(1) = vp(1)+M(5);
vp(2) = vp(2)+M(6);
vp(3) = vp(3)+M(7);
vp(4) = vp(4)+M(8);
vs(1) = vs(1)+M(9);
vs(2) = vs(2)+M(10);
vs(3) = vs(3)+M(11);
vs(4) = vs(4)+M(12);

































