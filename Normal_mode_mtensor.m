% ���������еĵ�����������
% �ӵѿ�������ϵת����������ϵ�����Ƕ�����Դ������ϵ��ָ��֮���ڱ仯��
% �ѿ�������ϵ�Ķ��壺x�ķ����Ǳ���y�ķ����Ƕ���z�ķ���������
% ������ϵ�Ķ��壺r�Ǿ���t����γ�ȷ���p�Ƕ����ȷ�����ʱ��
%
function [Msph] = Normal_mode_mtensor(Mcart)
Mxx = Mcart(1,1); %Mcart�ǵѿ���ϵ���������
Myy = Mcart(2,2);
Mzz = Mcart(3,3);
Mxy = Mcart(1,2);
Mxz = Mcart(1,3);
Myz = Mcart(2,3);
Mrr = Mzz;
Mrt = Mxz;
Mrp = -Myz;
Mtt = Mxx;
Mpp = Myy;
Mtp = -Mxy;
Msph = [Mrr Mrt Mrp;Mrt Mtt Mtp;Mrp Mtp Mpp]; %Msph��������ϵ���������
end