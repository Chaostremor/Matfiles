% 读取sac文件，并根据其中头段变量的值筛选台站，留下可用的台站数据
% 比如这里的筛选条件是：震中距dist要在 [30,90]度之间
%
clear ; clc ;
cd('C:\Users\Song Chao\Documents\MATLAB\myself\MatSAC') ;
datadir = 'G:\newzld\newzld_hinet\' ;   % 数据所在目录
%mkdir('G:\newzld\newzld_hinet_sort\') ;
sortdir = 'G:\newzld\newzld_hinet_sort\' ;   % 筛选后的目录
fid = fopen(strcat(datadir,'filelist')) ;      % strcat用于字符串连接，filelist事先编辑好
filelist = textscan(fid, '%s \n') ;  % 得到的是一个cell array,调用语句为flielist{1}  
%
fclose(fid) ;
flist = char(filelist{1}) ;   % 文件名字符串数组  
[sa,sb] = size(flist) ;    % 文件名字符串数组的大小
evlat = 42.7568 ;   % 事件纬度
evlon = 173.0766 ;   % 事件经度
%
j=1 ;
stlon = zeros(sa,1) ;
stlat =  zeros(sa,1) ;
disdeg = zeros(sa,1) ;
azideg = zeros(sa,1) ;
for i = 1:sa
    filename = strcat(datadir,flist(i,:)) ;  % 获取文件路径
    [t, data, SAChdr] = fget_sac(filename) ;   % 读取sac文件，t是时间，data是数据，SAChdr是头段变量
    stlon(i) = SAChdr.station.stlo ;   % 注意头段变量的存储方式
    stlat(i) = SAChdr.station.stla ;
    [disdeg(i), azideg(i)] = distance (evlat, evlon, stlat(i), stlon(i)) ;  % 两种算法有微小的差异
    %[diskm, disdeg, bazideg, azideg] = distaz (stlat, stlon, evlat, evlon) ;
    % 判断条件
    if (disdeg(i) >= 30.0) && (disdeg(i) <= 90.0)
        fremlist(j,:) = flist(i,:) ;
        disremin(j,1) = disdeg(i) ;
        copyfile(filename, strcat(sortdir, fremlist(j,:))) ;   % 复制文件到筛选后的目录
        j=j+1;
    end
end