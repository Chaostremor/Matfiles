% 将大震的所有子震源解的结果叠加在一起
%
clear; close all;
%
% 读入所有的结果文件，以竖直方向上的位移为例
uz1 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p1.uz');
uz2 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p2.uz');
uz3 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p3.uz');
uz4 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p4.uz');
uz5 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p5.uz');
uz6 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p6.uz');
uz7 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p7.uz');
uz8 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p8.uz');
uz9 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p9.uz');
uz10 = importdata('G:\qssp2010-code+input\server_result_Japan\Japan_total\Japan_out_p10.uz');
% 读入接收台站的所有信息
fid = fopen('G:\qssp2010-code+input\server_result_Japan\Japan_total\118stainfo.txt');
stainfo = textscan(fid, '%f %f %s %d %f ');
fclose(fid);
% 得到的是一个cell array,调用语句stainfo{1}、{2} ... , 以及stainfo{1}（1）。  
%
% 把波形数据和时间数据分开，只留下波形,并得到T_origin
Twin = uz1.data(end,1) - uz1.data(1,1);
dt = uz1.data(2,1) - uz1.data(1,1);
wavedatauz(1,:,:) = uz1.data(:,2:end);
wavedatauz(2,:,:) = uz2.data(:,2:end);
wavedatauz(3,:,:) = uz3.data(:,2:end);
wavedatauz(4,:,:) = uz4.data(:,2:end);
wavedatauz(5,:,:) = uz5.data(:,2:end);
wavedatauz(6,:,:) = uz6.data(:,2:end);
wavedatauz(7,:,:) = uz7.data(:,2:end);
wavedatauz(8,:,:) = uz8.data(:,2:end);
wavedatauz(9,:,:) = uz9.data(:,2:end);
wavedatauz(10,:,:) = uz10.data(:,2:end);
%wavedatauz是一个三维数组，第一维是文件个数，第二维是数据点数，第三维是台站个数
%wavedatauz(1,:,:) = uz1.data(:,2:end);
%
T_origin(1) = uz1.data(1,1);
T_origin(2) = uz2.data(1,1);
T_origin(3) = uz3.data(1,1);
T_origin(4) = uz4.data(1,1);
T_origin(5) = uz5.data(1,1);
T_origin(6) = uz6.data(1,1);
T_origin(7) = uz7.data(1,1);
T_origin(8) = uz8.data(1,1);
T_origin(9) = uz9.data(1,1);
T_origin(10) = uz10.data(1,1);
T_origin_max = max(T_origin);
% T_origin是一维数组，存放每个文件的起始时刻
%
[mz,nz,pz] = size(wavedatauz);
Treduction = stainfo{4};  % 开始计算的时刻
dist = stainfo{5};  % 台站震中距
%
% 为了把某一台站同一分量的所有文件相加，需要参照T_origin补零成同样size，对应相加
dim = ceil((T_origin_max + Twin) / dt) + 1;  
% dim是存放任意台站数据经补零之后的大小
uzi = zeros(mz,dim,pz);
% uzi是三维数组，第一维是文件数，第二维是dim，第三维是台站数，先置零，后面只用加上有数据的就可以
%
for k = 1:pz
for j = 1:mz
    %uzi(j,1:T_origin(j)/dt,k) = 0;
    uzi(j,ceil(T_origin(j)/dt)+1:ceil((Twin+T_origin(j))/dt)+1,k) = wavedatauz(j,:,k);
    % 以最小的T_origin为参考，其余的根据T_origin差值进行数组下标的位移
    %uzi(j,(Twin+T_origin(j))/dt+2:,k) = 0;
end
end
%uzik = uzi(:,:,1)';
%aa = sum(uzik,2);
uztotal = zeros(dim,pz);
for k = 1:pz
    uzik = uzi(:,:,k)';  % uzik用来临时存放某一台站的数据，转置使得数据每一列是一个文件的数据
    uztotal(:,k) = sum(uzik,2);  % uztotal是某一台站所有数据的求和，对uzik的列求和
end
figure (1)
for k = 1:pz
    %lengthk = double((Twin + Treduction(k)));
    %tuzk = linspace(0,lengthk,ceil(lengthk/dt)+1);
    tuzk = (double(Treduction(k)):dt:double(Twin+T_origin_max+Treduction(k))+dt)';
    % 冒号表达式构造时间向量，从Treduction开始，长度为Twin+T_origin_max
    uzkplot = uztotal(:,k)./(max(abs(uztotal(:,k)))) + double(dist(k))/2;
    % 用振幅绝对值的最大值作归一化，纵轴上以震中距大小区分
    plot(tuzk,uzkplot);hold on;
end
set(gca,'XTick',0:50:1200);
set(gca,'XLim',[0 1200]);  % 设置x轴的间隔，长度
