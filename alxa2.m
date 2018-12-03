clear;
lateq = 39.8;
loneq = 106.3;
workdir = 'G:\Alxa\';
fid = fopen(strcat(workdir,'X2Chninfo-675.txt'));
stainfo = textscan(fid,'%s %s %f %f %s %s %s %s %s %s %s %s %s %s %s %s');
nknm = char(stainfo{1});
stnm = char(stainfo{2});
lonst = stainfo{3};
latst = stainfo{4};
% ni = size(latst);
[dist, azi] = distance(lateq, loneq, latst, lonst) ;
dist = 111.0 *dist ;
[ni , mi] = size(dist);
j = 1;
for i = 1:ni
    if dist(i) >=50.0 && dist(i) <=400.0
        slecstnm(j,:) = stnm(i,:);
        slecdist(j) = dist(i);
        slecazi(j) = azi(i);
        j = j+1;
    end
end
[ni , mi] = size(slecdist);
for i = 5:10:355
    k = 1;
    for j = 1:ni
        if abs(slecazi(j)-i) <=1.0
            srt1azi(k) = slecazi(j);
            srt1stnm(k,:) = slecstnm(a,:);
            srt1dist(k) = slecazi(a);
    end
end

    
    
[srt1dist, ind1] = sort(slecdist);
srt1dist = srt1dist';
[ni , mi] = size(srt1dist);
for i=1:ni
    a= ind1(i);
    srt1stnm(i,:) = slecstnm(a,:);
    srt1azi(i) = slecazi(a);
end
srt1azi = srt1azi';
%[srt2azi, ind2] = sort(srt1azi);
[ni , mi] = size(srt1azi);
k=1;
%sadasd= find (abs(srt1azi-60) == min(abs(srt1azi-60)));
for j = 5:10:355
        index(k) = find (abs(srt1azi-j) == min(abs(srt1azi-j)));
        srt2azi(k) = srt1azi(index(k));
        srt2stnm(k,:) = srt1stnm(index(k),:);
        srt2dist(k) = srt1dist(index(k));
        k=k+1;
end
index = index';
%srtindex = sort(index);
srt2azi = srt2azi';
srt2dist = srt2dist';
% for i = 1:k-2
%     for j = i:k-1
%         if index(i) == index(j)
%             srt2azi(i)=[];
%             srt2stnm(i,:)=[];
%             srt2dist(i)=[];
%         end
%     end
% end
%[ni , mi] = size (srt2azi);
for i = 1:k-2
    if abs(srt2azi(i)-srt2azi(i+1)) <5.0
         srt2azi(i)=[];
         srt2stnm(i,:)=[];
         srt2dist(i)=[];
    end
end
[finaldist, ind2] = sort(srt2dist);
[ni , mi] = size(finaldist);
for i=1:ni
    a= ind2(i);
    finalstnm(i,:) = srt2stnm(a,:);
    finalazi(i) = srt2azi(a);
end
finalazi = finalazi';
