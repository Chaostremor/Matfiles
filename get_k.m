function [kkk,datam]=get_k(name,tend)
%name='/Users/Zhipeng/Documents/MUSIC/nepal/aileen_data/eu/eu_grid105stations10s0.5HzTo2Hz/HFdots_tcro'
data=load(name);
pp=find(data(:,1)<tend);
datainterp3=min(data(pp,3)):0.1:max(data(pp,3));
datainterp2=interp1(data(pp,3),data(pp,2),datainterp3);
%p=polyfit(datainterp3,datainterp2,1)
p=polyfix(datainterp3,datainterp2,1,[84.7314],[28.2305]);
%p=polyfix(data(pp,3),data(pp,2),1,[84.7314],[28.2305])
datam=data(pp,:);
kkk=p(1)


end
