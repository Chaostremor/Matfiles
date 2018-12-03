% to mexico 3d low-freq. result
% get bp rupture direction, speed, length and so on, plot then
% 
% Author: 
%     C. Song, 2017.10.4
%

%% initial setting
clear;
close all;
% initial parameter, change when needed
%%%%%%%%%%%%%%%%%%%
lon0=-93.715;           % lat0, lon0 denote epicenter
lat0=15.0678;
tend=60;                     % ending time
data = load('G:\BackProjection\mexico\AL0.2-1\mexico_al2183dstations10s0.2HzTo1Hz\HFdots_tc');          % HFdots dir
%%%%%%%%%%%%%%%%%%%

%%
t = data(:, 1);                                 % time, 1st col.
ind=find(t<tend);
t = t(ind);
lat = data(ind, 2);                          % latitude, 2nd col.
lon = data(ind, 3);                         % longitude, 2nd col.
power = data(ind, 4);                    % normalized power, 0-1
x = deg2km(lat - lat0);                          % coordinate transfer
y = deg2km((lon - lon0)*cosd(lat0));               % x-->north-->lat, y-->east-->lon
    
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y, x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on

% das = sqrt(x.^2+y.^2);                 % to test whether two algorithms is valid
% das2 = deg2km(distance(lat0, lon0, lat, lon));

% from the plot and summary_tc fig. , it is obvious that rupture only has
% one stage, with one uniform direction
% 
%%
% 
yinter = min(y): 1: max(y);         % loninter, latinter used to polyfix a line 
xinter = interp1(y, x, yinter);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(yinter, xinter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on
coef = polyfix(yinter, xinter, 1, 0, 0);            % polyfit points and must through 0,0, get slope
plot(yinter, polyval(coef, yinter), 'b-');
axis equal;
aaa = -20;                                % any value denote lon at the propagation side 
bbb = polyval(coef, aaa);         % y=polyval(slope, x) 
angle = atan2d(bbb, aaa);         % atan2d(Y,X), start from x-axis, clockwise 0~-180, counter 0~180
direction = 90-angle;               % rupture direction
if direction<0
    direction = direction+360;
end
%
yseq = min(y): 0.1: 0;             % lat1, lon1 used to plot
xseq = polyval(coef, yseq);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y , x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq, xseq, 'b-', 'linewidth', 2);
axis equal;
%
ang = direction - 90;               % rotate angle
R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];             % coordinate rotate matrix
yrot = zeros((length(y)), 1);
xrot = zeros((length(y)), 1);
dist = zeros((length(y)), 1);
for i = 1: length(y)
    new = R*[y(i); x(i)];          % coordinate transform
    yrot(i) = new(1, 1);
    xrot(i) = new(2, 1);
    if yrot(i)>0
        dist(i) = yrot(i);           % dist1 = positive coord.
    end
end
% dist = abs(y3);
figure
plot(t, dist, 'ko', 'MarkerFaceColor','r', 'markersize', 6);
%
ndata = [t x y dist power lat lon];
ndata = sortrows(ndata, 1);           % the times are calibrated, it may be integer and in ascending order 
t = ndata(:, 1);
x = ndata(:, 2);
y = ndata(:, 3);
dist = ndata(:, 4);
power = ndata(:, 5);
lat = ndata(:, 6);
lon = ndata(:, 7);

% no limit, keep all points
scale = max(dist);          % rupture length
tinter = min(t): 0.1: max(t);
distint = interp1(t, dist, tinter);
% coef1 = polyfix(tinter, distint, 1, t(1), dist(1));
coef = polyfit(tinter, distint, 1);
figure
plot(t, dist, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter, polyval(coef, tinter), 'b-');
speed = coef(1);                % rupture speed

% constraint with hypothesis: no retreat on one side
% this is for 210 stations
% distuse = [dist(1)];
% xuse = [x(1)];
% yuse = [y(1)];
% tuse = [t(1)];
% latuse = [lat(1)];
% lonuse = [lon(1)];
% powuse = [power(1)];
% for i = 2: length(y)
%     if dist(i)-max(dist(1: i-1)) >= -10
%         distuse = [distuse; dist(i)];
%         xuse = [xuse; x(i)];
%         yuse = [yuse; y(i)];
%         tuse = [tuse; t(i)];
%         latuse = [latuse; lat(i)];
%         lonuse = [lonuse; lon(i)];
%         powuse = [powuse; power(i)];
%     end
% end

% this is for 218 stations
distcoef = polyval(coef, t);
index=[];
for i = 1: length(y)
    if (dist(i)-distcoef(i) >= -10) && (dist(i)-distcoef(i) <= 20)
        index = [index; i];
    end
end
distuse = dist(index);
xuse = x(index);
yuse = y(index);
tuse = t(index);
latuse = lat(index);
lonuse = lon(index);
powuse = power(index);

scale1 = max(distuse);          % rupture length
tinter1 = min(tuse): 0.1: max(tuse);
distint1 = interp1(tuse, distuse, tinter1);
% coef2 = polyfix(tinter, distint, 1, 0, 0);
coef = polyfit(tinter1, distint1, 1);
figure
plot(tuse, distuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter1, polyval(coef, tinter1), 'b-');
speed1 = coef(1);                % rupture speed


%%
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(yuse, xuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(yseq, xseq, 'b-', 'linewidth', 2);
axis equal;

fdata = [tuse latuse lonuse powuse];
save('G:\BackProjection\mexico\AL0.2-1\mexico_al2183dstations10s0.2HzTo1Hz\HFdots_final_al3d_lf', 'fdata','-ascii');
save('G:\BackProjection\mexico\AL\mexico_al210stations10s0.5HzTo2Hz\al3dlf.mat');