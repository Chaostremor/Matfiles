% ���������������Դ��Ľ��������һ��
%
clear; close all;
%
% �������еĽ���ļ�������ֱ�����ϵ�λ��Ϊ��
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
% �������̨վ��������Ϣ
fid = fopen('G:\qssp2010-code+input\server_result_Japan\Japan_total\118stainfo.txt');
stainfo = textscan(fid, '%f %f %s %d %f ');
fclose(fid);
% �õ�����һ��cell array,�������stainfo{1}��{2} ... , �Լ�stainfo{1}��1����  
%
% �Ѳ������ݺ�ʱ�����ݷֿ���ֻ���²���,���õ�T_origin
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
%wavedatauz��һ����ά���飬��һά���ļ��������ڶ�ά�����ݵ���������ά��̨վ����
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
% T_origin��һά���飬���ÿ���ļ�����ʼʱ��
%
[mz,nz,pz] = size(wavedatauz);
Treduction = stainfo{4};  % ��ʼ�����ʱ��
dist = stainfo{5};  % ̨վ���о�
%
% Ϊ�˰�ĳһ̨վͬһ�����������ļ���ӣ���Ҫ����T_origin�����ͬ��size����Ӧ���
dim = ceil((T_origin_max + Twin) / dt) + 1;  
% dim�Ǵ������̨վ���ݾ�����֮��Ĵ�С
uzi = zeros(mz,dim,pz);
% uzi����ά���飬��һά���ļ������ڶ�ά��dim������ά��̨վ���������㣬����ֻ�ü��������ݵľͿ���
%
for k = 1:pz
for j = 1:mz
    %uzi(j,1:T_origin(j)/dt,k) = 0;
    uzi(j,ceil(T_origin(j)/dt)+1:ceil((Twin+T_origin(j))/dt)+1,k) = wavedatauz(j,:,k);
    % ����С��T_originΪ�ο�������ĸ���T_origin��ֵ���������±��λ��
    %uzi(j,(Twin+T_origin(j))/dt+2:,k) = 0;
end
end
%uzik = uzi(:,:,1)';
%aa = sum(uzik,2);
uztotal = zeros(dim,pz);
for k = 1:pz
    uzik = uzi(:,:,k)';  % uzik������ʱ���ĳһ̨վ�����ݣ�ת��ʹ������ÿһ����һ���ļ�������
    uztotal(:,k) = sum(uzik,2);  % uztotal��ĳһ̨վ�������ݵ���ͣ���uzik�������
end
figure (1)
for k = 1:pz
    %lengthk = double((Twin + Treduction(k)));
    %tuzk = linspace(0,lengthk,ceil(lengthk/dt)+1);
    tuzk = (double(Treduction(k)):dt:double(Twin+T_origin_max+Treduction(k))+dt)';
    % ð�ű��ʽ����ʱ����������Treduction��ʼ������ΪTwin+T_origin_max
    uzkplot = uztotal(:,k)./(max(abs(uztotal(:,k)))) + double(dist(k))/2;
    % ���������ֵ�����ֵ����һ���������������о��С����
    plot(tuzk,uzkplot);hold on;
end
set(gca,'XTick',0:50:1200);
set(gca,'XLim',[0 1200]);  % ����x��ļ��������
