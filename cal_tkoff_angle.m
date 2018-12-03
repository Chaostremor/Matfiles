% ������Դ��ih
%
% r0�ǵ���뾶
% h����Դ���ڵ���ȣ�rh=r0-h
% v0�ǵ��洦�ĵ��𲨲��٣����ݲ������ࣩ
% vh����Դ���ĵ������
% i0�������
% ��Գƽ��ʵ�Snell�����ǣ�rh*sin(ih)/vh=r0*sin(i0)/v0
% ����ǿ��ɵ��𲨵�ʱ������t��delta��h����ã�
% sin(i0)=v0*dt/ddelta������delta�Գ�������kmΪ��λ
% �����տ��Եõ���
% sin(ih)=r0/(r0-h)*vh*dt/ddelta
function [ih] = cal_tkoff_angle(r0,h,vh,dtddelta)
sinih = r0/(r0-h)*vh*dtddelta;
ih = asin(sinih);
ih = rad2deg(ih);
return
end