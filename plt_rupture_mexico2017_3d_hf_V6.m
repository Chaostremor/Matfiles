% to mexico 3d high-freq. result, version 6.0
% get bp rupture direction, speed, length and so on, plot then
% 
% Author: 
%     C. Song, 2018.3.30
% Revised based on v5
%    2018.4.8
%    2018.4.20
%    corresponds to the new results from runteleBPv3 or 3dtimelibv3 (different methods to determine the global maximum)
%    specificly for AL-neic location

%% initial setting
clear;
close all;
% initial parameter, change when needed
%%%%%%%%%%%%%%%%%%%
lon0 = -93.9067;           % lat0, lon0 denote epicenter
lat0 = 15.0356;
data = load('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\HFdots_tc');          % HFdots dir
%%%%%%%%%%%%%%%%%%%

%% determine the end of stacking time.
tstack = data(:, end);
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
x = deg2km(lat - lat0);                          % coordinate transfer
y = deg2km((lon - lon0)*cosd(lat0));               % x-->north-->lat, y-->east-->lon
    
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y, x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
axis equal;

figure
plot(trup, power, 'linestyle', '-', 'color', 'k', 'linewidth', 2);
%coef = polyfix(y, x, 1, 0, 0);
%plot(y, polyval(coef, y), 'b-');

% das = sqrt(x.^2+y.^2);                 % to test whether two algorithms is valid
% das2 = deg2km(distance(lat0, lon0, lat, lon));

% from the plot and summary_tc fig. , it is obvious that rupture can be
% divided into two stages, solve the direction and speed respectively

%%
% stage1, 1-32
corind1 = 32;
t1 =  [trup(1: 20); trup(22: 24); trup(26: corind1)];
lat1 = [lat(1: 20); lat(22: 24); lat(26: corind1)];
lon1 = [lon(1: 20); lon(22: 24); lon(26: corind1)];
power1 = [power(1: 20); power(22: 24); power(26: corind1)];
alpha1 = [alpha(1: 20); alpha(22: 24); alpha(26: corind1)];
x1 = [x(1: 20); x(22: 24); x(26: corind1)];
y1 = [y(1: 20); y(22: 24); y(26: corind1)];
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
coef11 = polyfix(x1inter, y1inter, 1, 0, 0);
plot(polyval(coef11, x1inter), x1inter, 'b-');
axis equal;
angle=atand(coef11(1));
if coef11(1)>=0       % rupture to right, look from the lower plane
    direction1 = angle;
else                           %  rupture to left
    direction1 = 360+angle;
end
xseq1 = 0: 0.1: max(x1);             % lat1, lon1 used to plot
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
% %
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
for i = 1: length(y1)
    new = R*[y1(i); x1(i)];          % coordinate transform
    y1rot(i) = new(1, 1);
    x1rot(i) = new(2, 1);
    if y1rot(i)>0
        dist1(i) = y1rot(i);           % dist1 = positive coord.
    end
end
% dist = abs(y3);
figure
plot(t1, dist1, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
%
ndata = [t1 x1 y1 dist1 power1 lat1 lon1 alpha1];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t1 = ndata(:, 1);
x1 = ndata(:, 2);
y1 = ndata(:, 3);
dist1 = ndata(:, 4);
power1 = ndata(:, 5);
lat1 = ndata(:, 6);
lon1 = ndata(:, 7);
alpha1 = ndata(:, 8);

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
corind2 = 35;
end2 = 46;
t2 = trup(corind2: end2);
lat2 = lat(corind2: end2);
lon2 = lon(corind2: end2);
power2 = power(corind2: end2);
x2 = x(corind2: end2);
y2 = y(corind2: end2);
alpha2 = alpha(corind2: end2);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y2, x2, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(y2(1), x2(1),'ko', 'MarkerFaceColor','b', 'markersize', 6); hold on
plot(y2(end), x2(end), 'ko', 'MarkerFaceColor','y', 'markersize', 6); hold on
axis equal;
%

[y2sort, index]=sort(y2,'ascend');
x2sort=x2(index);
y2inter = (y2sort(1): 0.1: y2sort(end))';       % loninter, latinter used to polyfix a line
x2inter = interp1(y2, x2, y2inter);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y2inter, x2inter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on

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
xseq2 = min(x2): 0.1: max(x2);             % lat1, lon1 used to plot
yseq2 = polyval(coef21, xseq2);
%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%
% coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));            % polyfit points and must through 0,0, get slope
% % coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));
% plot(y2inter, polyval(coef21, y2inter), 'b-');
% axis equal;
% angle=atand(coef21(1));
% if y2(end)>=y2(1)       % rupture to right, look from the lower plane
%     direction2= 90-angle;
% else                           %  rupture to left
%     direction2 = 270-angle;
% end
% yseq2 = min(y2): 0.1: y2(1);             % lat1, lon1 used to plot
% xseq2 = polyval(coef21, yseq2);
% %%%%%%%%%%%%%%%%%%%%%%%

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
for i = 1: length(y2)
    new = R*[y2(i)-y2(1); x2(i)-x2(1)];          % coordinate transform
    y2rot(i) = new(1, 1);
    x2rot(i) = new(2, 1);
    if y2rot(i)>0
        dist2(i) = y2rot(i);           % dist1 = positive coord.
    end
end
% dist = abs(y3);
figure
plot(t2, dist2, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
%
ndata = [t2 x2 y2 dist2 power2 lat2 lon2 alpha2];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t2 = ndata(:, 1);
x2 = ndata(:, 2);
y2 = ndata(:, 3);
dist2 = ndata(:, 4);
power2 = ndata(:, 5);
lat2 = ndata(:, 6);
lon2 = ndata(:, 7);
alpha2 = ndata(:, 8);

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
distuse = [distuse1; distuse2];

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(yuse, xuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq1, xseq1, 'b-', 'linewidth', 2);
plot(yseq2, xseq2, 'b-', 'linewidth', 2);
axis equal;

figure
plot(tuse, powuse, 'linestyle', '-', 'color', 'k', 'linewidth', 2);

fdata = [tuse latuse lonuse powuse alphause distuse xuse yuse];
save('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\HFdots_final_al3d_hf', 'fdata','-ascii');
save('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\al3dhf.mat');

%% convert cartisan coordinate to map projection
xpro1 = min(x1)-4: 2: max(x1);             % lat1, lon1 used to plot
ypro1 = polyval(coef11, xpro1);
xpro2 = min(x2): 2: max(x2)+2;             % lat1, lon1 used to plot
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
save('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\profile1', 'profile1','-ascii');
save('G:\BackProjection\mexico\AL-neic\mexico_al_neic2103dstations10s0.5HzTo2Hz\profile2', 'profile2','-ascii');

%%
dep0 = 56.67;
dep = 70;
dip = 83;
disu = dep0/tand(dip);
x2stu = xpro2(1)-disu*sind(360-direction2);        % x-->north-->lat, y-->east-->lon
y2stu = ypro2(1)-disu*cosd(360-direction2);
x2edu = xpro2(end)-disu*sind(360-direction2);
y2edu = ypro2(end)-disu*cosd(360-direction2);

lat2stu = km2deg(x2stu)+lat0;
lon2stu = km2deg(y2stu/cosd(lat0))+lon0;
lat2edu = km2deg(x2edu)+lat0;
lon2edu = km2deg(y2edu/cosd(lat0))+lon0;

disd = (dep-dep0)/tand(dip);
x2std = xpro2(1)+disd*sind(360-direction2);        % x-->north-->lat, y-->east-->lon
y2std = ypro2(1)+disd*cosd(360-direction2);
x2edd = xpro2(end)+disd*sind(360-direction2);
y2edd = ypro2(end)+disd*cosd(360-direction2);

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
plot([y2std; y2edd], [x2std; x2edd], 'r-', 'linewidth', 2); hold on
plot(ypro2, xpro2, 'b-', 'linewidth', 2);
axis equal

