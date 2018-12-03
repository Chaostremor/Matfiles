% get the crustal info of study region, (102-109E, 36-42N)
% epicenter location is 39.7663, 106.4304

clear;
[lon, lat] = meshgrid(102: 1: 109, 36: 1: 42);
c = getcrust(lat, lon);
[nlat, nlon] = size(lat);
% tksed = reshape(-c.top(:, 6), [nlat nlon]);
tksed = reshape(c.thk(:, 3), [nlat nlon])+reshape(c.thk(:, 4), [nlat nlon])+reshape(c.thk(:, 5), [nlat nlon]);
tkcrtu = reshape(c.thk(:, 6), [nlat nlon]);
tkcrtm = reshape(c.thk(:, 7), [nlat nlon]);
tkcrtl = reshape(c.thk(:, 8), [nlat nlon]);
tkmoho = reshape(-c.top(:, 9), [nlat nlon]);
avetksed = mean(tksed(:));
avetkcrtu = mean(tkcrtu(:));
avetkcrtm = mean(tkcrtm(:));
avetkcrtl = mean(tkcrtl(:));
avetkmoho = mean(tkmoho(:));

vpsed = reshape(c.vp(:, 5), [nlat nlon]);
vpcrtu = reshape(c.vp(:, 6), [nlat nlon]);
vpcrtm = reshape(c.vp(:, 7), [nlat nlon]);
vpcrtl = reshape(c.vp(:, 8), [nlat nlon]);
vpmoho = reshape(c.vp(:, 9), [nlat nlon]);
avevpsed = mean(vpsed(:));
avevpcrtu = mean(vpcrtu(:));
avevpcrtm = mean(vpcrtm(:));
avevpcrtl = mean(vpcrtl(:));
avevpmoho = mean(vpmoho(:));

vssed = reshape(c.vs(:, 5), [nlat nlon]);
vscrtu = reshape(c.vs(:, 6), [nlat nlon]);
vscrtm = reshape(c.vs(:, 7), [nlat nlon]);
vscrtl = reshape(c.vs(:, 8), [nlat nlon]);
vsmoho = reshape(c.vs(:, 9), [nlat nlon]);
avevssed = mean(vssed(:));
avevscrtu = mean(vscrtu(:));
avevscrtm = mean(vscrtm(:));
avevscrtl = mean(vscrtl(:));
avevsmoho = mean(vsmoho(:));

avemod(1, 1) = avetksed;
avemod(2, 1) = avetkcrtu;
avemod(3, 1) = avetkcrtm;
avemod(4, 1) = avetkcrtl;
avemod(5, 1) = avetkmoho;
avemod(1, 2) = avevpsed;
avemod(2, 2) = avevpcrtu;
avemod(3, 2) = avevpcrtm;
avemod(4, 2) = avevpcrtl;
avemod(5, 2) = avevpmoho;
avemod(1, 3) = avevssed;
avemod(2, 3) = avevscrtu;
avemod(3, 3) = avevscrtm;
avemod(4, 3) = avevscrtl;
avemod(5, 3) = avevsmoho;

%%

% lon0 = 106.4304;
% lat0 = 39.7663;
lon0 = 108;
lat0 = 37.2;
b = getcrust(lat0, lon0);

% vpwat0 = b.vp(1);
% vpice0 = b.vp(2);
% vpsedu0 = b.vp(3);
% vpsedm0 = b.vp(4);
% vpsedl0 = b.vp(5);
% vswat0 = b.vs(1);
% vsice0 = b.vs(2);
% vssedu0 = b.vs(3);
% vssedm0 = b.vs(4);
% vssedl0 = b.vs(5);

% epi(1,1) = b.thk(1);
% epi(2,1) = b.thk(2);
% epi(3,1) = b.thk(3);
% epi(4,1) = b.thk(4);
% epi(5,1) = b.thk(5);
% epi(6,1) = b.thk(6);
% epi(7,1) = b.thk(7);
% epi(8,1) = b.thk(8);
% epi(9,1) = -b.top(9);
% epi(1, 2) = b.vp(1);
% epi(2, 2) = b.vp(2);
% epi(3, 2) = b.vp(3);
% epi(4, 2) = b.vp(4);
% epi(5, 2) = b.vp(5);
% epi(6, 2) = b.vp(6);
% epi(7, 2) = b.vp(7);
% epi(8, 2) = b.vp(8);
% epi(9, 2) = b.vp(9);
% epi(1, 3) = b.vs(1);
% epi(2, 3) = b.vs(2);
% epi(3, 3) = b.vs(3);
% epi(4, 3) = b.vs(4);
% epi(5, 3) = b.vs(5);
% epi(6, 3) = b.vs(6);
% epi(7, 3) = b.vs(7);
% epi(8, 3) = b.vs(8);
% epi(9, 3) = b.vs(9);
epi(1: 9, 1) = cat(2, b.thk(1: 8), -b.top(9));
epi(1: 9, 2) = b.vp(1: 9);
epi(1: 9, 3) = b.vs(1: 9);








