% �����������������������ò���������λ��ͶӰ��P��ƫ�����ϣ�
% ���ϳɵ����Ĵ���������
% ��̨վ�ķ�λ����phi��ֱ��P�����������inc
% ÿһ����������un��ue��uz���ϳɵ�u�����ǵĹ�ϵ�ǣ�
% u=un*cos(phi)*sin(inc)+ue*sin(phi)*sin(inc)-uz*cos(inc)
%
clear; close all;
%
% ����mrr������λ��
mrruz = importdata('G:\qssp2010-code+input\Azimuth9\mrr\azi_out_mrr.uz');
mrrue = importdata('G:\qssp2010-code+input\Azimuth9\mrr\azi_out_mrr.ue');
mrrun = importdata('G:\qssp2010-code+input\Azimuth9\mrr\azi_out_mrr.un');
% ����mtt������λ��
mttuz = importdata('G:\qssp2010-code+input\Azimuth9\mtt\azi_out_mtt.uz');
mttue = importdata('G:\qssp2010-code+input\Azimuth9\mtt\azi_out_mtt.ue');
mttun = importdata('G:\qssp2010-code+input\Azimuth9\mtt\azi_out_mtt.un');
% ����mpp������λ��
mppuz = importdata('G:\qssp2010-code+input\Azimuth9\mpp\azi_out_mpp.uz');
mppue = importdata('G:\qssp2010-code+input\Azimuth9\mpp\azi_out_mpp.ue');
mppun = importdata('G:\qssp2010-code+input\Azimuth9\mpp\azi_out_mpp.un');
% ����mrt������λ��
mrtuz = importdata('G:\qssp2010-code+input\Azimuth9\mrt\azi_out_mrt.uz');
mrtue = importdata('G:\qssp2010-code+input\Azimuth9\mrt\azi_out_mrt.ue');
mrtun = importdata('G:\qssp2010-code+input\Azimuth9\mrt\azi_out_mrt.un');
% ����mrp������λ��
mrpuz = importdata('G:\qssp2010-code+input\Azimuth9\mrp\azi_out_mrp.uz');
mrpue = importdata('G:\qssp2010-code+input\Azimuth9\mrp\azi_out_mrp.ue');
mrpun = importdata('G:\qssp2010-code+input\Azimuth9\mrp\azi_out_mrp.un');
% ����mtp������λ��
mtpuz = importdata('G:\qssp2010-code+input\Azimuth9\mtp\azi_out_mtp.uz');
mtpue = importdata('G:\qssp2010-code+input\Azimuth9\mtp\azi_out_mtp.ue');
mtpun = importdata('G:\qssp2010-code+input\Azimuth9\mtp\azi_out_mtp.un');
% ����ʱ������
Tun = mrruz.data(:,1);
% z�������������
uz(1,:,:) = mrruz.data(:,2:end);
uz(2,:,:) = mttuz.data(:,2:end);
uz(3,:,:) = mppuz.data(:,2:end);
uz(4,:,:) = mrtuz.data(:,2:end);
uz(5,:,:) = mrpuz.data(:,2:end);
uz(6,:,:) = mtpuz.data(:,2:end);
% e�������������
ue(1,:,:) = mrrue.data(:,2:end);
ue(2,:,:) = mttue.data(:,2:end);
ue(3,:,:) = mppue.data(:,2:end);
ue(4,:,:) = mrtue.data(:,2:end);
ue(5,:,:) = mrpue.data(:,2:end);
ue(6,:,:) = mtpue.data(:,2:end);
% n�������������
un(1,:,:) = mrrun.data(:,2:end);
un(2,:,:) = mttun.data(:,2:end);
un(3,:,:) = mppun.data(:,2:end);
un(4,:,:) = mrtun.data(:,2:end);
un(5,:,:) = mrpun.data(:,2:end);
un(6,:,:) = mtpun.data(:,2:end);
%
[ms,ns,ps] = size(uz);
% 
x = 0;   % x=phi
y = 20.01;   % y=inc
cosx = cos(deg2rad(x));
sinx = sin(deg2rad(x));
cosy = cos(deg2rad(y));
siny = sin(deg2rad(y));
%
u = zeros(size(uz));
for i=1:ms
    for j=1:ps
        for k=1:ns
            % u=un*cos(phi)*sin(inc)+ue*sin(phi)*sin(inc)-uz*cos(inc)
            u(i,k,j) = un(i,k,j)*cosx*siny + ue(i,k,j)*sinx*siny -uz(i,k,j)*cosy;
        end
    end
end























