% coded by Z. P. Liu, modified by C. Song 

clear
close all
lat=28.2305;        % epicenter lat and lon
lon=84.7314;
tend=55;             % end time
xextend=1.4        % lon range
% HFdots format
% time(calibrated) lat lon power
[kau,dataau]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/au/AU_grid983dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend    % discretize lon
y=kau*(x-lon)+lat;
plot(x,y,'b','LineWidth',2)     % fig after interpolating HF point

hold on
plot(dataau(:,3),dataau(:,2),'bo','MarkerFaceColor','b','MarkerSize',2)    % data ---> actual HF point
for i=1:length(dataau)
   ddd=distance(dataau(i,2),dataau(i,3),y,x);           % distance between points, (lat1, lon1, lat2, lon2) 
   index=find(ddd==min(ddd));
   distau(i)=deg2km(distance(lat,lon,y(index),x(index)))   % distau is distance from epicenter to HF point
end

[keu,dataeu]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/eu/eu_grid1053dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend
y=keu*(x-lon)+lat;
plot(x,y,'g','LineWidth',2)
hold on
plot(dataeu(:,3),dataeu(:,2),'go','MarkerFaceColor','g','MarkerSize',2)

for i=1:length(dataeu)
   ddd=distance(dataeu(i,2),dataeu(i,3),y,x);
   index=find(ddd==min(ddd));
   disteu(i)=deg2km(distance(lat,lon,y(index),x(index)))
end

[kus,dataus]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/alaska/AL_grid1163dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend
y=kus*(x-lon)+lat;
plot(x,y,'r','LineWidth',2)
hold on
plot(dataus(:,3),dataus(:,2),'ro','MarkerFaceColor','r','MarkerSize',2)
set(gca,'DataAspectRatio',[1/cosd(lat) 1 1])                       % set axis x, y ratio,   x, y, z
plot(lon,lat,'Kp','MarkerSize',10,'MarkerFaceColor','k')
for i=1:length(dataus)
   ddd=distance(dataus(i,2),dataus(i,3),y,x);
   index=find(ddd==min(ddd));
   distus(i)=deg2km(distance(lat,lon,y(index),x(index)))
end
ylim([27.7,28.4])
xlim([84.5 ,86.2])
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
 chsize(20)
 
 figure
 plot(dataau(:,1),distau,'bo','MarkerFaceColor','b')
 
 hold on
 plot(dataeu(:,1),disteu,'go','MarkerFaceColor','g')
 plot(dataus(:,1),distus,'ro','MarkerFaceColor','r')
  
 vau=polyfit(dataau(:,1),distau',1)
 veu=polyfit(dataeu(:,1),disteu',1)
 vus=polyfit(dataus(:,1),distus',1)
 
 
 ll=15;
 arrow([40,20],[40+ll,20+vau*ll],'Color','b','Width',1)
 hold on
 arrow([40,20],[40+ll,20+veu*ll],'Color','g','Width',1)
 hold on
 arrow([40,20],[40+ll,20+vus*ll],'Color','r','Width',1)