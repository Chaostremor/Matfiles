% to mexico 3d high-freq. result, version 7.0
% get bp rupture direction, speed, length and so on, plot then
% 
% Author: 
%     C. Song, 2018.3.30
% Revised based on v5
%    2018.4.9
%    corresponds to the new results from runteleBPv3 or 3dtimelibv3 (different methods to determine the global maximum)
%    specificly for AL-usgs location

%% initial setting
clear;
close all;
% initial parameter, change when needed
%%%%%%%%%%%%%%%%%%%
lon0 = -93.899;           % lat0, lon0 denote epicenter
lat0 = 15.022;
data = load('G:\BackProjection\mexico\AL-usgs\mexico_al_usgs2103dstations10s0.5HzTo2Hz\HFdots_tc');          % HFdots dir
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
corind1=32;
t1 =  [trup(1: 19); trup(21: 22); trup(24); trup(26: corind1)];
lat1 = [lat(1: 19); lat(21: 22); lat(24); lat(26: corind1)];
lon1 = [lon(1: 19); lon(21: 22); lon(24); lon(26: corind1)];
power1 = [power(1: 19); power(21: 22); power(24); power(26: corind1)];
alpha1 = [alpha(1: 19); alpha(21: 22); alpha(24); alpha(26: corind1)];
x1 = [x(1: 19); x(21: 22); x(24); x(26: corind1)];
y1 = [y(1: 19); y(21: 22); y(24); y(26: corind1)];

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y1, x1, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
% 

[y1sort, index]=sort(y1,'ascend');
x1sort=x1(index);
y1inter = (y1sort(1): 0.1: y1sort(end))';       % loninter, latinter used to polyfix a line
x1inter = interp1(y1, x1, y1inter);

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y1inter, x1inter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on
coef11 = polyfix(y1inter, x1inter, 1, 0, 0);            % polyfit points and must through 0,0, get slope
plot(y1inter, polyval(coef11, y1inter), 'b-');
axis equal;
angle=atand(coef11(1));
if y1(end)>=y1(1)       % rupture to right, look from the lower plane
    direction1 = 90-angle;
else                           %  rupture to left
    direction1 = 270-angle;
end
%
yseq1 = min(y1): 0.1: 0;             % lat1, lon1 used to plot
xseq1 = polyval(coef11, yseq1);
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
ndata = [t1 x1 y1 dist1 power1 lat1 lon1];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t1 = ndata(:, 1);
x1 = ndata(:, 2);
y1 = ndata(:, 3);
dist1 = ndata(:, 4);
power1 = ndata(:, 5);
lat1 = ndata(:, 6);
lon1 = ndata(:, 7);

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
for i = 2: length(y1)
    if dist1(i)-max(dist1(1: i-1)) >= -15
        distuse1 = [distuse1; dist1(i)];
        xuse1 = [xuse1; x1(i)];
        yuse1 = [yuse1; y1(i)];
        tuse1 = [tuse1; t1(i)];
        latuse1 = [latuse1; lat1(i)];
        lonuse1 = [lonuse1; lon1(i)];
        powuse1 = [powuse1; power1(i)];
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
corind2 = 34;
end2 = 46;
t2 =  [trup(corind2: end2)];
lat2 = [lat(corind2: end2)];
lon2 = [lon(corind2: end2)];
power2 = [power(corind2: end2)];
x2 = [x(corind2: end2)];
y2 = [y(corind2: end2)];
% t2 =  [trup(corind2); trup(35: end2)];
% lat2 = [lat(corind2); lat(35: end2)];
% lon2 = [lon(corind2); lon(35: end2)];
% power2 = [power(corind2); power(35: end2)];
% x2 = [x(corind2); x(35: end2)];
% y2 = [y(corind2); y(35: end2)];
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
coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));            % polyfit points and must through 0,0, get slope
% coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));
plot(y2inter, polyval(coef21, y2inter), 'b-');
axis equal;

% y2inter = min(y2): 0.2: max(y2);
% % y2inter = min(y2): 0.2: max(y2);         % loninter, latinter used to polyfix a line 
% x2inter = interp1(y2, x2, y2inter);
% figure
% plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
% plot(y2inter, x2inter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on
% coef21 = polyfix(y2inter, x2inter, 1, y2(1), x2(1));            % polyfit points and must through 0,0, get slope
% plot(y2inter, polyval(coef21, y2inter), 'b-');
% axis equal;
angle=atand(coef21(1));
if y2(end)>=y2(1)       % rupture to right, look from the lower plane
    direction2= 90-angle;
else                           %  rupture to left
    direction2 = 270-angle;
end

%
yseq2 = min(y2): 0.1: y2(1);             % lat1, lon1 used to plot
xseq2 = polyval(coef21, yseq2);
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
ndata = [t2 x2 y2 dist2 power2 lat2 lon2];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t2 = ndata(:, 1);
x2 = ndata(:, 2);
y2 = ndata(:, 3);
dist2 = ndata(:, 4);
power2 = ndata(:, 5);
lat2 = ndata(:, 6);
lon2 = ndata(:, 7);

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
for i = 2: length(y2)
     if dist2(i)-max(dist2(1: i-1)) >= -15
%      if dist2(i)-max(dist2(1: i-1)) >= -50
        distuse2 = [distuse2; dist2(i)];
        xuse2 = [xuse2; x2(i)];
        yuse2 = [yuse2; y2(i)];
        tuse2 = [tuse2; t2(i)];
        latuse2 = [latuse2; lat2(i)];
        lonuse2 = [lonuse2; lon2(i)];
        powuse2 = [powuse2; power2(i)];
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

figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(yuse, xuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq1, xseq1, 'b-', 'linewidth', 2);
plot(yseq2, xseq2, 'b-', 'linewidth', 2);
axis equal;

figure
plot(tuse, powuse, 'linestyle', '-', 'color', 'k', 'linewidth', 2);

fdata = [tuse latuse lonuse powuse];
save('G:\BackProjection\mexico\AL-usgs\mexico_al_usgs2103dstations10s0.5HzTo2Hz\HFdots_final_al3d_hf', 'fdata','-ascii');
save('G:\BackProjection\mexico\AL-usgs\mexico_al_usgs2103dstations10s0.5HzTo2Hz\al3dhf.mat');

figure
for i=1: length(tuse)
    marker='o';
    scatter(lonuse(i), latuse(i), powuse(i)*100, marker); hold on
end
axis equal;
plot(lon0,lat0,'r*','MarkerSize',10); hold on
set(gca,'DataAspectRatio',[1/cosd(lat0) 1 1]);
ux=linspace(-1+lat0,2.5+lat0,35);
uy=linspace(-2+lon0,2+lon0,40);
ylim([min(ux) max(ux)]);
xlim([min(uy) max(uy)]);