% coded by Z. P. Liu, modified by C. Song 

clear
close all
lat=28.2305;           % epicenter lat and lon
lon=84.7314;
tend=55;                % end time
xextend=1.4;          % lon range
% HFdots format
% time(calibrated) lat lon power
[kau,dataau]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/au/AU_grid98stations10s0.5HzTo2Hz/HFdots_tc',tend);

x=lon:0.01:lon+xextend;            % discretize lon
y=kau*(x-lon)+lat;
subplot(2,2,1)

pau=plot(x,y,'m','LineWidth',2)       % fig after interpolating HF point, pau is handle
hold on
plot(dataau(:,3),dataau(:,2),'mo','MarkerFaceColor','m','MarkerSize',2)          % data ---> actual HF point
for i=1:length(dataau)
   ddd=distance(dataau(i,2),dataau(i,3),y,x);             % distance between points, (lat1, lon1, lat2, lon2) 
   index=find(ddd==min(ddd));
   distau(i)=deg2km(distance(lat,lon,y(index),x(index)));       % distau is distance from epicenter to HF point
end

[keu,dataeu]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/eu/eu_grid105stations10s0.5HzTo2Hz/HFdots_tc',tend);

x=lon:0.01:lon+xextend;
y=keu*(x-lon)+lat;
peu=plot(x,y,'b','LineWidth',2)
hold on
plot(dataeu(:,3),dataeu(:,2),'bo','MarkerFaceColor','b','MarkerSize',2)
for i=1:length(dataeu)
   ddd=distance(dataeu(i,2),dataeu(i,3),y,x);
   index=find(ddd==min(ddd));
   disteu(i)=deg2km(distance(lat,lon,y(index),x(index)));
end


[kus,dataus]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/alaska/AL_grid116stations10s0.5HzTo2Hz/HFdots_tc',tend);

x=lon:0.01:lon+xextend;
y=kus*(x-lon)+lat;
pus=plot(x,y,'g','LineWidth',2)
hold on
plot(dataus(:,3),dataus(:,2),'go','MarkerFaceColor','g','MarkerSize',2)
set(gca,'DataAspectRatio',[1/cosd(lat) 1 1])
plot(lon,lat,'rp','MarkerSize',15,'MarkerFaceColor','r')         % plot epicenter 
for i=1:length(dataus)
   ddd=distance(dataus(i,2),dataus(i,3),y,x);
   index=find(ddd==min(ddd));
   distus(i)=deg2km(distance(lat,lon,y(index),x(index)));
end
ylim([27.7,28.4])
xlim([84.5 ,86.2])
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
%title('(a)')
  legend([peu,pus,pau],'EU','NA','AU','location','southwest')
 chsize(15)
  subplot(2,2,3)

 
 plot(dataau(:,1),distau,'mo','MarkerFaceColor','m','MarkerSize',2)
 
 hold on
 plot(dataeu(:,1),disteu,'bo','MarkerFaceColor','b','MarkerSize',2)
 plot(dataus(:,1),distus,'go','MarkerFaceColor','g','MarkerSize',2)
 vau=polyfit(dataau(:,1),distau',1)
 veu=polyfit(dataeu(:,1),disteu',1)
 vus=polyfit(dataus(:,1),distus',1)
  ll=15;
 arrow([40,20],[40+ll,20+vau*ll],8,'Color','m','Width',1)
 hold on
 arrow([40,20],[40+ll,20+veu*ll],8,'Color','b','Width',1)
 hold on
 arrow([40,20],[40+ll,20+vus*ll],8,'Color','g','Width',1)
 
 %%
 clear

lat=28.2305;
lon=84.7314;
tend=55;
xextend=1.4;
[kau,dataau]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/au/AU_grid983dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend;
y=kau*(x-lon)+lat;
xlabel('Time (s)');
ylabel('Distance (km)');

 % title('(c)')
   chsize(15)
subplot(2,2,2)

plot(x,y,'m','LineWidth',2)

hold on
plot(dataau(:,3),dataau(:,2),'mo','MarkerFaceColor','m','MarkerSize',2)
for i=1:length(dataau)
   ddd=distance(dataau(i,2),dataau(i,3),y,x);
   index=find(ddd==min(ddd));
   distau(i)=deg2km(distance(lat,lon,y(index),x(index)));
end

[keu,dataeu]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/eu/eu_grid1053dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend;
y=keu*(x-lon)+lat;
plot(x,y,'b','LineWidth',2)
hold on
plot(dataeu(:,3),dataeu(:,2),'bo','MarkerFaceColor','b','MarkerSize',2)

for i=1:length(dataeu)
   ddd=distance(dataeu(i,2),dataeu(i,3),y,x);
   index=find(ddd==min(ddd));
   disteu(i)=deg2km(distance(lat,lon,y(index),x(index)));
end

[kus,dataus]=get_k('/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/alaska/AL_grid1163dstations10s0.5HzTo2Hz/HFdots_tcro',tend);

x=lon:0.001:lon+xextend;
y=kus*(x-lon)+lat;
plot(x,y,'g','LineWidth',2)
hold on
plot(dataus(:,3),dataus(:,2),'go','MarkerFaceColor','g','MarkerSize',2)
set(gca,'DataAspectRatio',[1/cosd(lat) 1 1])
plot(lon,lat,'rp','MarkerSize',15,'MarkerFaceColor','r')
for i=1:length(dataus)
   ddd=distance(dataus(i,2),dataus(i,3),y,x);
   index=find(ddd==min(ddd));
   distus(i)=deg2km(distance(lat,lon,y(index),x(index)));
end
ylim([27.7,28.4])
xlim([84.5 ,86.2])
xlabel('Longitude (^o)');
ylabel('Latitude (^o)');
 %title('(b)')
 chsize(15)

subplot(2,2,4)

 plot(dataau(:,1),distau,'mo','MarkerFaceColor','m','MarkerSize',2)
 
 hold on
 plot(dataeu(:,1),disteu,'bo','MarkerFaceColor','b','MarkerSize',2)
 plot(dataus(:,1),distus,'go','MarkerFaceColor','g','MarkerSize',2)
  
 vau=polyfit(dataau(:,1),distau',1)
 veu=polyfit(dataeu(:,1),disteu',1)
 vus=polyfit(dataus(:,1),distus',1)
 
 
 ll=15;
pau=arrow([40,20],[40+ll,20+vau*ll],8,'Color','m','Width',1)
 hold on
 peu=arrow([40,20],[40+ll,20+veu*ll],8,'Color','b','Width',1)
 hold on
 pus=arrow([40,20],[40+ll,20+vus*ll],8,'Color','g','Width',1)
 xlabel('Time (s)');
ylabel('Distance (km)');

%title('(d)')
  chsize(15)