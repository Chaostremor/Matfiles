% compare results before and after relocation
%
% C. Song, 
%
% Last modified: 2018/01/23.

clear; close all;

% load data after relocation
% load tlag between real and syn and real pg/sg arrival time
% load 

weidir = 'G:\Alxa\nodecimate\3test400\' ;    % 数据所在目录
fid1 = fopen(strcat(weidir,'nweight.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnm = char(weight{1});     % 台站名是第一列
odist = weight{2};
[sa, sb] = size(stnm);
maxdist = 200;
for ii=1:sa
    if odist(ii) > maxdist        
        index = ii-1;
        break
    end
end

load('alxa47_18_pg_tlag.mat');
oplaguse = cat(1, totallag(2: 11), totallag(13: index));
synpguse = cat(1, synpg(2: 11), synpg(13: index));
realpguse = cat(1, realpg(2: 11), realpg(13: index));
odistuse = cat(1, odist(2: 11), odist(13: index));
stnmuse = cat(1, stnm(2: 11, :), stnm(13: index, :));

orelaplag = oplaguse/synpguse;

load('alxa47_18_sg_tlag.mat');
oslaguse = cat(1, totallag(2: 11), totallag(13: index));
realsguse = cat(1, realsg(2: 11), realsg(13: index));

weidir = 'G:\Alxa\alxa49_18\' ;    % 数据所在目录
fid1 = fopen(strcat(weidir,'reweight_reloc.dat')) ;      % strcat用于字符串连接
weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid1) ;
stnmcheck = char(weight{1});
ndist = weight{2};

figure (1)
plot (odistuse, oplaguse, 'b.', 'Markersize', 20);
set(gca, 'YLIM', [-4, 4]);
set(gca, 'YTick', -4: 2: 4);

% dorigin = sum(oplaguse)/length(oplaguse);
% rmeanlag = oplaguse - dorigin;
% orelaplag = rmeanlag/synpguse;

figure (2)
plot (odistuse, orelaplag, 'b.', 'Markersize', 20);
set(gca, 'YLIM', [-0.1, 0.1]);
set(gca, 'YTick', -0.1: 0.05: 0.1);

% figure (2)
% plot (odistuse, oslaguse, 'r.', 'Markersize', 20);
% set(gca, 'YLIM', [-4, 4]);
% set(gca, 'YTick', -4: 2: 4);


% model para.  -->  alxa47_18
vs = [3.28; 3.55; 4.18; 4.48];
vp = [5.50; 5.96; 6.60; 8.12];
tk= [6.7; 31.4; 15.6; 0.0];
dep = 18;

% syn Pg time 
nsample1 = 10000;
rayp1 = linspace(0, 1.0/vp(2)-0.0000005, nsample1)';
xpg = zeros(nsample1,1);
tpg = zeros(nsample1,1);
for ii = 1: nsample1
    etap1 = sqrt((1.0/vp(1))^2-rayp1(ii)^2);
    etap2 = sqrt((1.0/vp(2))^2-rayp1(ii)^2);
    xpg(ii) = rayp1(ii) * ( tk(1)/etap1 + (dep-tk(1))/etap2 );
    tpg(ii) = rayp1(ii)*xpg(ii) + tk(1)*etap1 + (dep-tk(1))*etap2;                              
end
synpg = interp1(xpg, tpg, ndist, 'spline');       % model determined pg curve
nplag = realpguse - synpg;
% nrelaplag = nplag/synpg;

% syn Sg time
nsample1 = 10000;
rays1 = linspace(0, 1.0/vs(2)-0.0000005, nsample1)';
xsg = zeros(nsample1,1);
tsg = zeros(nsample1,1);
for ii = 1: nsample1
    etas1 = sqrt((1.0/vs(1))^2-rays1(ii)^2);
    etas2 = sqrt((1.0/vs(2))^2-rays1(ii)^2);
    xsg(ii) = rays1(ii) * ( tk(1)/etas1 + (dep-tk(1))/etas2 );
    tsg(ii) = rays1(ii)*xsg(ii) + tk(1)*etas1 + (dep-tk(1))*etas2;                              
end
synsg = interp1(xsg, tsg, ndist, 'spline');
nslag = realsguse - synsg;

figure (3)
plot (ndist, nplag, 'b.', 'Markersize', 20);
set(gca, 'YLIM', [-4, 4]);
set(gca, 'YTick', -4: 2: 4);


% figure (4)
% plot (ndist, nslag, 'b.', 'Markersize', 20);
% set(gca, 'YLIM', [-4, 4]);
% set(gca, 'YTick', -4: 2: 4);


dorigin = sum(nplag)/length(nplag);
rmeanlag = nplag - dorigin;
nrelaplag = rmeanlag/synpg;

figure (4)
plot (ndist, nrelaplag, 'b.', 'Markersize', 20);
set(gca, 'YLIM', [-0.1, 0.1]);
set(gca, 'YTick', -0.1: 0.05: 0.1);

oorigin = 27.210;
norigin = oorigin +dorigin;

save('relocation_check.mat'); 