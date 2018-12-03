% img=imread('1.png');
% figure(1)
% imshow(img);
% hold on
% 
% [x0,y0]=ginput(2);
% 
% for i=1: 30
%     [x(i),y(i)]=ginput(1);
%     plot(x(i),y(i),'ob','MarkerFaceColor','b');
% end
    
% lblon = 104.5;
% lblat = 38;
% rtlon = 108.5;
% rtlat = 41;
% dlat = (rtlat-lblat)/(y0(1)-y0(2));
% dlon = (rtlon-lblon)/(x0(2)-x0(1));
% lon = lblon+(x-x0(1))*dlon;
% lat = rtlat-(y-y0(2))*dlat;
% 
% f1lon = lon(1: 10);
% f1lat = lat(1: 10);
% f2lon = lon(11: 20);
% f2lat = lat(11: 20);
% f3lon = lon(21: end);
% f3lat = lat(21: end);
% figure
% plot(f1lon, f1lat); hold on
% plot(f2lon, f2lat); hold on
% plot(f3lon, f3lat); hold on
% 
% f1 = [f1lon' f1lat'];
% f2 = [f2lon' f2lat'];
% f3 = [f3lon' f3lat'];

img=imread('3.png');
figure
imshow(img);
hold on

[x0,y0]=ginput(2);

npt = 15;
for i=1: npt
    [x1(i), y1(i)]=ginput(1);
    plot(x1(i),y1(i),'ob','MarkerFaceColor','b');
end
lblon = 165;
lblat = 50;
rtlon = 175;
rtlat = 55;
dlat = (rtlat-lblat)/(y0(1)-y0(2));
dlon = (rtlon-lblon)/(x0(2)-x0(1));
lon = lblon+(x1-x0(1))*dlon;
lat = rtlat-(y1-y0(2))*dlat;

f4lon = lon(1: npt);
f4lat = lat(1: npt);
figure
plot(f4lon, f4lat);
f4 = [f4lon' f4lat'];




