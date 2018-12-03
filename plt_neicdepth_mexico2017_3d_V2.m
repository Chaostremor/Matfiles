% to mexico 3d high-freq. result with depth, version 2.0
% get bp rupture direction, speed, length and so on, plot then
% 
% Author: 
%     C. Song, 2018.5.12
% Revised based on v1
% 
%    corresponds to the new results from runteleBPv3 or 3dtimelibv3 (different methods to determine the global maximum)
%    specificly for AL-neic location with depth variation
%
%    first, shift the first bp dot (and others) to the epicenter; second, estimate the rupture direction, fix at the epicenter. 
%    OR, estimate the rupture direction, NO fix at the epicenter. 
%    Both way should obtain the same rupture direction

%% initial setting
clear;
close all;
% initial parameter, change when needed
%%%%%%%%%%%%%%%%%%%
lon0 = -93.9067;           % lat0, lon0 denote epicenter
lat0 = 15.0356;
data = load('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\HFdots_tc');          % HFdots dir
%%%%%%%%%%%%%%%%%%%

%% determine the end of stacking time.
tstack = data(:, 6);
power = data(:, 4);                    % normalized power, 0-1
coeff =data(:, 5);
alpha=power.*coeff;
figure
plot(tstack, alpha, 'linestyle', '-', 'color', 'k', 'linewidth', 2); hold on
line([0 80], [0.1 0.1], 'linewidth', 0.5, 'color', [180/255 180/255 180/255], 'linestyle', '--'); hold on

%%
ind=1:47;         % ending index, determined from alpha figure
trup = data(ind, 1);                                 % time, 1st col.
lat = data(ind, 2);                          % latitude, 2nd col.
lon = data(ind, 3);                         % longitude, 2nd col.
power = data(ind, 4);                    % normalized power, 0-1
coeff =data(ind, 5);                       % cross-correlation coefficients, 0-1
bdep = data(ind, 7);
x = deg2km(lat - lat0);                          % coordinate transfer
y = deg2km((lon - lon0)*cosd(lat0));               % x-->north-->lat, y-->east-->lon
    
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y, x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
axis equal;

figure
plot(trup, power, 'linestyle', '-', 'color', 'k', 'linewidth', 2);

% shift all points relative to the epicenter
lat = lat-(lat(1)-lat0);
lon = lon-(lon(1)-lon0);
x = deg2km(lat - lat0);                          % coordinate transfer
y = deg2km((lon - lon0)*cosd(lat0));               % x-->north-->lat, y-->east-->lon
    
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y, x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
axis equal;

% from the plot and summary_tc fig. , it is obvious that rupture can be
% divided into two stages, solve the direction and speed respectively

%%
% stage1, 1-32
corind1 = 32;
% t1 =  [trup(1: 5); trup(6: corind1)];
% lat1 = [lat(1: 5); lat(6: corind1)];
% lon1 = [lon(1: 5); lon(6: corind1)];
% power1 = [power(1: 5); power(6: corind1)];
% alpha1 = [alpha(1: 5); alpha(6: corind1)];
% x1 = [x(1: 5); x(6: corind1)];
% y1 = [y(1: 5); y(6: corind1)];
% bdep1= [bdep(1: 5); bdep(6: corind1)];

t1 =  [trup(1: 6); trup(9: 10); trup(12: 25); trup(27: 28); trup(30: corind1)];
lat1 = [lat(1: 6); lat(9: 10); lat(12: 25); lat(27: 28); lat(30: corind1)];
lon1 = [lon(1: 6); lon(9: 10); lon(12: 25); lon(27: 28); lon(30: corind1)];
power1 = [power(1: 6); power(9: 10); power(12: 25); power(27: 28); power(30: corind1)];
alpha1 = [alpha(1: 6); alpha(9: 10); alpha(12: 25); alpha(27: 28); alpha(30: corind1)];
x1 = [x(1: 6); x(9: 10); x(12: 25); x(27: 28); x(30: corind1)];
y1 = [y(1: 6); y(9: 10); y(12: 25); y(27: 28); y(30: corind1)];
bdep1=[bdep(1: 6); bdep(9: 10); bdep(12: 25); bdep(27: 28); bdep(30: corind1)];
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y1, x1, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(y1(1), x1(1), 'ko', 'MarkerFaceColor','b', 'markersize', 6); hold on
plot(y1(end), x1(end), 'ko', 'MarkerFaceColor','y', 'markersize', 6); hold on
% 
% [y1sort, index]=sort(y1,'ascend');
% x1sort=x1(index);
% y1inter = (y1sort(1): 0.1: y1sort(end))';       % loninter, latinter used to polyfix a line
% x1inter = interp1(y1, x1, y1inter);

[x1sort, index]=sort(x1,'ascend');
y1sort=y1(index);
x1inter = (x1sort(1): 0.1: x1sort(end))';         % loninter, latinter used to polyfix a line 
y1inter = interp1(x1sort, y1sort, x1inter);

% y1inter = min(y1): 0.1: max(y1);         % loninter, latinter used to polyfix a line 
% x1inter = interp1(y1, x1, y1inter);

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y1inter, x1inter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on

%%%%%%%%%%%%%%%%%%%%
%%% polyfix a line along latitude (x)
coef11 = polyfix(x1inter, y1inter, 1, 0, 0);          % fix at the epicenter
% coef11 = polyfix(x1inter, y1inter, 1, x1(1), y1(1));   % fix at the start point of stage 1 
% coef11 = polyfit(x1inter, y1inter, 1);                  % no fix at one point
plot(polyval(coef11, x1inter), x1inter, 'b-');
axis equal;
angle=atand(coef11(1));
if coef11(1)>=0       % rupture to right, look from the lower plane
    direction1 = angle;
else                           %  rupture to left
    direction1 = 360+angle;
end
xseq1 = 0: 1: max(x1);             % lat1, lon1 used to plot
yseq1 = polyval(coef11, xseq1);
%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%
% %%% polyfix a line along longitude (y)
% coef11 = polyfix(y1inter, x1inter, 1, 0, 0);            % polyfit points and must through 0,0, get slope
% plot(y1inter, polyval(coef11, y1inter), 'b-');
% axis equal;
% angle=atand(coef11(1));
% if y1(end)>=y1(1)       % rupture to right, look from the lower plane
%     direction1 = 90-angle;
% else                           %  rupture to left
%     direction1 = 270-angle;
% end
% yseq1 = min(y1): 0.1: 0;             % lat1, lon1 used to plot
% xseq1 = polyval(coef11, yseq1);
% %%%%%%%%%%%%%%%%%%%%

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y1 , x1, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq1, xseq1, 'b-', 'linewidth', 2);
axis equal;
%
ang = direction1 - 90;               % rotate angle
R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];             % coordinate rotate matrix
y1rot = zeros((length(y1)), 1);
x1rot = zeros((length(y1)), 1);
dist1 = zeros((length(y1)), 1);
odist1 = zeros((length(y1)), 1);
for i = 1: length(y1)
    new = R*[y1(i); x1(i)];          % coordinate transform
    y1rot(i) = new(1, 1);
    x1rot(i) = new(2, 1);
    if y1rot(i)>0
        dist1(i) = y1rot(i);           % dist1 = positive coord.
    end
    odist1(i) = y1rot(i);
end
% dist = abs(y3);
figure
plot(t1, dist1, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
figure
plot(t1, odist1, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
%
ndata = [t1 x1 y1 dist1 power1 lat1 lon1 alpha1 bdep1 odist1];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t1 = ndata(:, 1);
x1 = ndata(:, 2);
y1 = ndata(:, 3);
dist1 = ndata(:, 4);
power1 = ndata(:, 5);
lat1 = ndata(:, 6);
lon1 = ndata(:, 7);
alpha1 = ndata(:, 8);
bdep1 = ndata(:, 9);
odist1 = ndata(:, 10);

% no limit, keep all points
scale1 = max(dist1);          % rupture length
tinter1 = min(t1): 0.1: max(t1);
distint1 = interp1(t1, dist1, tinter1);
% coef1 = polyfix(tinter, distint, 1, t(1), dist(1));
coef12 = polyfit(tinter1, distint1, 1);
figure
plot(t1, dist1, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter1, polyval(coef12, tinter1), 'b-');
speed12 = coef12(1);                % rupture speed

% constraint with hypothesis: no retreat on one side
distuse1 = [dist1(1)];
xuse1 = [x1(1)];
yuse1 = [y1(1)];
tuse1 = [t1(1)];
latuse1 = [lat1(1)];
lonuse1 = [lon1(1)];
powuse1 = [power1(1)];
alphause1 = [alpha1(1)];
bdepuse1 = [bdep1(1)];
odistuse1 = [odist1(1)];
for i = 2: length(y1)
    if dist1(i)-max(dist1(1: i-1)) >= -15
        distuse1 = [distuse1; dist1(i)];
        xuse1 = [xuse1; x1(i)];
        yuse1 = [yuse1; y1(i)];
        tuse1 = [tuse1; t1(i)];
        latuse1 = [latuse1; lat1(i)];
        lonuse1 = [lonuse1; lon1(i)];
        powuse1 = [powuse1; power1(i)];
        alphause1 = [alphause1; alpha1(i)];
        bdepuse1 = [bdepuse1; bdep1(i)];
        odistuse1 = [odistuse1; odist1(i)];
    end
end
scale11 = max(distuse1);          % rupture length
tinter11 = min(tuse1): 0.1: max(tuse1);
distint11 = interp1(tuse1, distuse1, tinter11);
% coef2 = polyfix(tinter, distint, 1, 0, 0);
coef13 = polyfit(tinter11, distint11, 1);
figure
plot(tuse1, distuse1, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter11, polyval(coef13, tinter11), 'b-');
speed13 = coef13(1);                % rupture speed

%%
% stage2, 36: end
corind2 = 33;
end2 = 45;
t2 =  [trup(corind2); trup(35: 37); trup(39: 41); trup(43: end2)];
lat2 = [lat(corind2); lat(35: 37); lat(39: 41); lat(43: end2)];
lon2 = [lon(corind2); lon(35: 37); lon(39: 41);lon(43: end2)];
power2 = [power(corind2); power(35: 37); power(39: 41); power(43: end2)];
alpha2 = [alpha(corind2); alpha(35: 37); alpha(39: 41); alpha(43: end2)];
x2 = [x(corind2); x(35: 37); x(39: 41); x(43: end2)];
y2 = [y(corind2); y(35: 37); y(39: 41); y(43: end2)];
bdep2 = [bdep(corind2); bdep(35: 37); bdep(39: 41); bdep(43: end2)];
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y2, x2, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(y2(1), x2(1),'ko', 'MarkerFaceColor','b', 'markersize', 6); hold on
plot(y2(end), x2(end), 'ko', 'MarkerFaceColor','y', 'markersize', 6); hold on
axis equal;
%

% [y2sort, index]=sort(y2,'ascend');
% x2sort=x2(index);
% y2inter = (y2sort(1): 0.1: y2sort(end))';       % loninter, latinter used to polyfix a line
% x2inter = interp1(y2, x2, y2inter);

[x2sort, index]=sort(x2,'ascend');
y2sort=y2(index);
x2inter = (x2sort(1): 1: x2sort(end))';         % loninter, latinter used to polyfix a line 
y2inter = interp1(x2sort, y2sort, x2inter);
 figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y2inter, x2inter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on
% coef21= polyfit(x2inter, y2inter, 1);
% plot(polyval(coef21, x2inter), x2inter, 'b-');

% coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));            % polyfit points and must through 0,0, get slope
% plot(y2inter, polyval(coef21, y2inter), 'b-');

%%%%%%%%%%%%%%%%%%%%%%%
coef21 = polyfix(x2inter, y2inter, 1, x2(1), y2(1));
plot(polyval(coef21, x2inter), x2inter, 'b-');
axis equal;
angle=atand(coef21(1));
if coef21(1)>=0       % rupture to right, look from the lower plane
    direction2 = angle;
else                           %  rupture to left
    direction2 = 360+angle;
end
xseq2 = min(x2): 1: max(x2);             % lat1, lon1 used to plot
yseq2 = polyval(coef21, xseq2);
%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y2 , x2, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq2, xseq2, 'b-', 'linewidth', 2);
axis equal;
%
ang = direction2 - 90;               % rotate angle
R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];             % coordinate rotate matrix
y2rot = zeros((length(y2)), 1);
x2rot = zeros((length(y2)), 1);
dist2 = zeros((length(y2)), 1);
odist2 = zeros((length(y2)), 1);
for i = 1: length(y2)
    new = R*[y2(i)-y2(1); x2(i)-x2(1)];          % coordinate transform
    y2rot(i) = new(1, 1);
    x2rot(i) = new(2, 1);
    if y2rot(i)>0
        dist2(i) = y2rot(i);           % dist1 = positive coord.
    end
    odist2(i) = y2rot(i);
end
% dist = abs(y3);
figure
plot(t2, dist2, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
figure
plot(t2, odist2, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
%
ndata = [t2 x2 y2 dist2 power2 lat2 lon2 alpha2 bdep2 odist2];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t2 = ndata(:, 1);
x2 = ndata(:, 2);
y2 = ndata(:, 3);
dist2 = ndata(:, 4);
power2 = ndata(:, 5);
lat2 = ndata(:, 6);
lon2 = ndata(:, 7);
alpha2 = ndata(:, 8);
bdep2 = ndata(:, 9);
odist2 = ndata(:, 10);

% no limit, keep all points
scale2 = max(dist2);          % rupture length
tinter2 = min(t2): 0.1: max(t2);
distint2 = interp1(t2, dist2, tinter2);
% coef1 = polyfix(tinter, distint, 1, t(1), dist(1));
coef22 = polyfit(tinter2, distint2, 1);
figure
plot(t2, dist2, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter2, polyval(coef22, tinter2), 'b-');
speed22 = coef22(1);                % rupture speed

% constraint with hypothesis: no retreat on one side
distuse2 = [dist2(1)];
xuse2 = [x2(1)];
yuse2 = [y2(1)];
tuse2 = [t2(1)];
latuse2 = [lat2(1)];
lonuse2 = [lon2(1)];
powuse2 = [power2(1)];
alphause2 = [alpha2(1)];
bdepuse2 = [bdep2(1)];
odistuse2 = [odist2(1)];
for i = 2: length(y2)
     if dist2(i)-max(dist2(1: i-1)) >= -15
%       if dist2(i)-max(dist2(1: i-1)) >= -50
        distuse2 = [distuse2; dist2(i)];
        xuse2 = [xuse2; x2(i)];
        yuse2 = [yuse2; y2(i)];
        tuse2 = [tuse2; t2(i)];
        latuse2 = [latuse2; lat2(i)];
        lonuse2 = [lonuse2; lon2(i)];
        powuse2 = [powuse2; power2(i)];
        alphause2 = [alphause2; alpha2(i)];
        bdepuse2 = [bdepuse2; bdep2(i)];
        odistuse2 = [odistuse2; odist2(i)];
     end
end

% % manually choose points
% distuse2 = [dist2(2: 3); dist2(5); dist2(7: end)];
% xuse2 = [x2(2: 3); x2(5); x2(7: end)];
% yuse2 = [y2(2: 3); y2(5); y2(7: end)];
% tuse2 = [t2(2: 3); t2(5); t2(7: end)];
% latuse2 = [lat2(2: 3); lat2(5); lat2(7: end)];
% lonuse2 = [lon2(2: 3); lon2(5); lon2(7: end)];
% powuse2 = [power2(2: 3); power2(5); power2(7: end)];

scale21 = max(distuse2);          % rupture length
tinter21 = min(tuse2): 0.1: max(tuse2);
distint21 = interp1(tuse2, distuse2, tinter21);
% coef2 = polyfix(tinter, distint, 1, 0, 0);
coef23 = polyfit(tinter21, distint21, 1);
figure
plot(tuse2, distuse2, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter21, polyval(coef23, tinter21), 'b-');
speed23 = coef23(1);                % rupture speed

%%
xuse = [xuse1; xuse2];
yuse = [yuse1; yuse2];
tuse = [tuse1; tuse2];
latuse = [latuse1; latuse2];
lonuse = [lonuse1; lonuse2];
powuse = [powuse1; powuse2];
alphause = [alphause1; alphause2];
bdepuse = [bdepuse1; bdepuse2];
distuse = [distuse1; distuse2];
odistuse2 = odistuse2+max(odistuse1);
odistuse = [odistuse1; odistuse2];

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(yuse, xuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq1, xseq1, 'b-', 'linewidth', 2);
plot(yseq2, xseq2, 'b-', 'linewidth', 2);
axis equal;

figure
plot(tuse, powuse, 'linestyle', '-', 'color', 'k', 'linewidth', 2);

fdata = [tuse latuse lonuse powuse alphause xuse yuse];
% depdata = [tuse powuse*0.7 bdepuse distuse odistuse];
depdata = [odistuse bdepuse tuse powuse distuse];
save('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\HFdots_al3d_depth_hf', 'fdata','-ascii');
save('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\HFdots_depth_profile', 'depdata','-ascii');
save('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\al3ddepthhf.mat');

%% convert cartisan coordinate to map projection
xpro1 = min(x1): 2: max(x1);             % lat1, lon1 used to plot
ypro1 = polyval(coef11, xpro1);
xpro2 = min(x2): 2: max(x2);             % lat1, lon1 used to plot
ypro2 = polyval(coef21, xpro2);

latpro1 = km2deg(xpro1)+lat0;
lonpro1 = km2deg(ypro1/cosd(lat0))+lon0;
latpro2 = km2deg(xpro2)+lat0;
lonpro2 = km2deg(ypro2/cosd(lat0))+lon0;

% x = deg2km(lat - lat0);                          % coordinate transfer
% y = deg2km((lon - lon0)*cosd(lat0));               % x-->north-->lat, y-->east-->lon
figure
plot(lon0, lat0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(lonpro1, latpro1, 'b-', 'linewidth', 2);
plot(lonpro2, latpro2, 'b-', 'linewidth', 2);
set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1]);
profile1 = [lonpro1' latpro1'];
profile2 = [lonpro2' latpro2'];
save('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\profile1', 'profile1','-ascii');
save('G:\BackProjection\mexico\depth\mexico_neic_depth2103dstations10s0.5HzTo2Hz\profile2', 'profile2','-ascii');

%%
dep0 = 56.67;
dep=70;
dip = 83;
disu = dep0/tand(dip);
x2stu = xpro2(1)+disu*sind(direction2);        % x-->north-->lat, y-->east-->lon
y2stu = ypro2(1)-disu*cosd(direction2);
x2edu = xpro2(end)+disu*sind(direction2);
y2edu = ypro2(end)-disu*cosd(direction2);

lat2stu = km2deg(x2stu)+lat0;
lon2stu = km2deg(y2stu/cosd(lat0))+lon0;
lat2edu = km2deg(x2edu)+lat0;
lon2edu = km2deg(y2edu/cosd(lat0))+lon0;

disd = (dep-dep0)/tand(dip);
x2std = xpro2(1)-disd*sind(direction2);        % x-->north-->lat, y-->east-->lon
y2std = ypro2(1)+disd*cosd(direction2);
x2edd = xpro2(end)-disd*sind(direction2);
y2edd = ypro2(end)+disd*cosd(direction2);

lat2std = km2deg(x2std)+lat0;
lon2std = km2deg(y2std/cosd(lat0))+lon0;
lat2edd = km2deg(x2edd)+lat0;
lon2edd = km2deg(y2edd/cosd(lat0))+lon0;

figure
plot(lon0, lat0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot([lon2stu; lon2edu], [lat2stu; lat2edu], 'g-', 'linewidth', 2);
plot([lon2std; lon2edd], [lat2std; lat2edd], 'r-', 'linewidth', 2);
plot(lonpro2, latpro2, 'b-', 'linewidth', 2);
set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1]);

figure
plot([y2stu; y2edu], [x2stu; x2edu], 'g-', 'linewidth', 2); hold on
plot([y2std; y2edd], [x2std; x2edd], 'r-', 'linewidth', 2);
plot(ypro2, xpro2, 'b-', 'linewidth', 2);
axis equal

fplane = [lon2stu lat2stu; lon2edu lat2edu; lonpro2(end) latpro2(end); lon2edd lat2edd; lon2std lat2std; lonpro2(1) latpro2(1); lon2stu lat2stu];
figure
plot(lon0, lat0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(fplane(:, 1), fplane(:, 2), 'g-', 'linewidth', 2);
set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1]);






