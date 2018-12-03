% get bp rupture direction, speed, length and so on, plot then
% Author: 
%     C. Song, 2017.7.28
%

%% initial setting
clear;
close all;
% initial parameter, change when needed
%%%%%%%%%%%%%%%%%%%
lon0=168.8148;           % lat0, lon0 denote epicenter
lat0=54.4715;
tend=70;                     % ending time
% data = load('G:\BackProjection\russia\EU\russia_eu263dstations10s0.5HzTo2Hz\HFdots_tc');
data = load('G:\BackProjection\russia\NA\russia_na2083dstations10s0.5HzTo2Hz\HFdots_tc');
%%%%%%%%%%%%%%%%%%%

%%
t = data(:, 1);
ind=find(t<tend);
t = t(ind);
lat = data(ind, 2);
lon = data(ind, 3);
power = data(ind, 4);
x = deg2km(lat - lat0);                          % coordinate transfer
y = deg2km((lon - lon0)*cosd(lat0));      % x-->north, y-->east
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y', 'markersize', 20); hold on
plot(y , x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
coef = polyfix(y, x, 1, 0, 0);
plot(y, polyval(coef, y), 'b-');

yinter = min(y): 1: max(y);         % loninter, latinter used to polyfix a line 
xinter = interp1(y, x, yinter);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(yinter, xinter, 'ko', 'MarkerFaceColor','r','markersize', 6); hold on
coef = polyfix(yinter, xinter, 1, 0, 0);
plot(yinter, polyval(coef, yinter), 'b-');
aaa = 120;
bbb = polyval(coef, aaa);
direction = atan2d(aaa, bbb);       % rupture direction

y1 = 0: 0.1: max(y);             % lat1, lon1 used to plot
x1 = polyval(coef, y1);
figure
plot(0, 0, 'kp', 'MarkerFaceColor','y','markersize', 20); hold on
plot(y , x, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(y1, x1, 'b-', 'linewidth', 2);
axis equal;

ang = direction - 90;
R = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];             % coordinate rotate matrix
y3 = zeros((length(y)), 1);
x3 = zeros((length(y)), 1);
dist = zeros((length(y)), 1);
for i = 1: length(y)
    new = R*[y(i); x(i)];          % coordinate transform
    y3(i) = new(1, 1);
    x3(i) = new(2, 1);
    if y3(i)>0
        dist(i) = y3(i);
    end
end
% dist = abs(y3);
figure
plot(t, dist, 'ko', 'MarkerFaceColor','r', 'markersize', 6);

ndata = [t x y dist power lat lon];
ndata = sortrows(ndata, 1);
t = ndata(:, 1);
x = ndata(:, 2);
y = ndata(:, 3);
dist = ndata(:, 4);
power = ndata(:, 5);
lat = ndata(:, 6);
lon = ndata(:, 7);

scale1 = max(dist);          % rupture length
tinter = min(t): 0.1: max(t);
distint = interp1(t, dist, tinter);
% coef1 = polyfix(tinter, distint, 1, t(1), dist(1));
coef1 = polyfit(tinter, distint, 1);
figure
plot(t, dist, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter, polyval(coef1, tinter), 'b-');
speed1 = coef1(1);                % rupture speed

%%
distuse = [];
xuse = [];
yuse = [];
tuse = [];
powuse = [];
latuse = [];
lonuse = [];
for i = 1: length(y)
%     distfit = (t(i)-t(1))*speed1+dist(1);
    distfit = polyval(coef1, t(i));
    if (dist(i)-distfit)>=-10 && (dist(i)-distfit)<=10
        distuse = [distuse; dist(i)];
        xuse = [xuse; x(i)];
        yuse = [yuse; y(i)];
        latuse = [latuse; lat(i)];
        lonuse = [lonuse; lon(i)];
        tuse = [tuse; t(i)];
        powuse = [powuse; power(i)];
    end
end
scale2 = max(distuse);          % rupture length
tinter = min(tuse): 0.1: max(tuse);
distint = interp1(tuse, distuse, tinter);
coef2 = polyfit(tinter, distint, 1);
figure
plot(tuse, distuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
plot(tinter, polyval(coef2, tinter), 'b-');
speed2 = coef2(1);                % rupture speed

%%
% % hypothesis: no retreat on one side
% distuse = [dist(1)];
% xuse = [x(1)];
% yuse = [y(1)];
% tuse = [t(1)];
% powuse = [power(1)];
% for i = 1: length(y)
%     if (dist(i)>= max(distuse(:))-5) && (dist(i)<= max(distuse(:))+10)
%         distuse = [distuse; dist(i)];
%         xuse = [xuse; x(i)];
%         yuse = [yuse; y(i)];
%         tuse = [tuse; t(i)];
%         powuse = [powuse; power(i)];
%     end
% end
% figure
% plot(tuse, distuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
% 
% scale2 = max(distuse);          % rupture length
% tinter = min(tuse): 0.1: max(tuse);
% distint = interp1(tuse, distuse, tinter);
% % coef2 = polyfix(tinter, distint, 1, 0, 0);
% coef2 = polyfit(tinter, distint, 1);
% figure
% plot(tuse, distuse, 'ko', 'MarkerFaceColor','r', 'markersize', 6); hold on
% plot(tinter, polyval(coef2, tinter), 'b-');
% speed2 = coef2(1);                % rupture speed

fdata = [tuse latuse lonuse powuse];
% save('G:\BackProjection\russia\EU\russia_eu263dstations10s0.5HzTo2Hz\HFdots_final_eu', 'fdata','-ascii');
save('G:\BackProjection\russia\NA\russia_na2083dstations10s0.5HzTo2Hz\HFdots_final_na', 'fdata','-ascii');