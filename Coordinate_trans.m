% ����ֱ������ϵ�º�������ϵ�µ��������������任
% Mxx,Myy,Mzz,Mxy,Mxz,Myz��ֱ������ϵ�µľ���������
% Mrr,Mtt,Mpp,Mrt,Mrp,Mtp��������ϵ�µľ���������
% �����x�Ǳ�����y�Ƕ�����z����ֱ����
% ͬʱ��r�Ǿ���t����γ��ָ��p�Ƕ���ָ����ʱ�룬��y������
% numda�Ǵ�y����ʱ��ľ��ȣ�0-360�㣩��phi����γ�ȣ�0-180�㣩
%
function [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp] = Coordinate_trans(Mxx,Myy,Mzz,Mxy,Mxz,Myz,numda,phi)
Mcart = [Mxx,Mxy,Mxz;Mxy,Myy,Myz;Mxz,Myz,Mzz]; %Mcart�ǵѿ���ϵ���������
a = deg2rad(numda);
b = deg2rad(phi);
cart2sph = [sin(a)*cos(b) cos(a)*cos(b) sin(b);-sin(a)*sin(b) -cos(a)*sin(b) cos(b);cos(a) -sin(a) 0];
%sph2cart = [sin(a)*cos(b) -sin(a)*sin(b) cos(a);cos(a)*cos(b) -cos(a)*sin(b) -sin(a);sin(b) cos(b) 0];
Msph = cart2sph * Mcart * (cart2sph)'; %Msph��������ϵ���������
%Mcart = sph2cart * Msph * (sph2cart)';
Mrr = Msph(1,1);
Mtt = Msph(2,2);
Mpp = Msph(3,3);
Mrt = Msph(1,2);
Mrp = Msph(1,3);
Mtp = Msph(2,3);
end