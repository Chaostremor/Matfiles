% 绘制某一个台站三个位移分量
%
clear; close all;
%
k = 1;
%
un = importdata('G:\qssp2010-code+input\Japantest\2Japan_out.un');
Tun = un.data(:,1);
%
Tun = Tun - 148 + 460;
%
wavedataun = un.data(:,2:end);
[mn,nn] = size(wavedataun);
unst1 = wavedataun(:,k);
%
ue = importdata('G:\qssp2010-code+input\Japantest\2Japan_out.ue');
Tue = ue.data(:,1);
%
Tue = Tue - 148 + 460;
%
wavedataue = ue.data(:,2:end);
[me,ne] = size(wavedataue);
uest1 = wavedataue(:,k);
%
uz = importdata('G:\qssp2010-code+input\Japantest\2Japan_out.uz');
Tuz = uz.data(:,1);
%
Tuz = Tuz - 148 + 460;
%
wavedatauz = uz.data(:,2:end);
[mz,nz] = size(wavedatauz);
uzst1 = wavedatauz(:,k);
%
treduction = load('G:\qssp2010-code+input\Japan\treduction.txt');
%
st1nm = un.colheaders{1,k+1};   % 读取struct的文件头型数据，注意用法
%
subplot(3,1,1);
plot(Tun,unst1,'-b');
xlabel('time /s');
ylabel('amplitude');
title(strcat('station name: ',st1nm,' , component: un'));
set(gca,'XTick',0:50:1500);
%
subplot(3,1,2);
plot(Tue,uest1,'-b');
xlabel('time /s');
ylabel('amplitude');
title(strcat('station name: ',st1nm,' , component: ue'));
set(gca,'XTick',0:50:1500);
%
subplot(3,1,3);
plot(Tuz,uzst1,'-b');
xlabel('time /s');
ylabel('amplitude');
title(strcat('station name: ',st1nm,' , component: uz'));
set(gca,'XTick',0:50:1500);
