%��֪һ���������еľ�γ�ȣ��Լ�̨վ�����еķ�λ�Ǻ����о�
%Ҫ��̨վ�ľ�γ��
%�����������ǵ����ҹ�ʽ����
%�����еľ�γ�ȷֱ�Ϊ��s_numda��s_phi������0-360��0-180��
%���о�Ϊdist(0-180)��
%��λ��Ϊazi(0-360)��
%�����̨վ�ľ�γ�ȷֱ�Ϊ��r_numda,r_phi������0-360��0-180��
%������õĸ��Ƕȵ�ȡֵ��Χ��ͬ����Ҫ���򵥵��޸�
%Last modified: 2016.10.15
%
function [r_numda,r_phi] = Receiver_loc_cal(s_numda,s_phi,dist,azi)
if azi<=180
   s_numda_rad = deg2rad(s_numda);
   s_phi_rad = deg2rad(s_phi);
   dist_rad = deg2rad(dist);
   azi_rad = deg2rad(azi);
   c_r_phi_rad = cos(s_phi_rad).*cos(dist_rad) + sin(s_phi_rad).*sin(dist_rad).*cos(azi_rad);
   r_phi_rad = acos(c_r_phi_rad);
   r_phi = rad2deg(r_phi_rad);
   %dnumda�Ǿ��Ȳ�,�㷨���ƣ�s_phi������0
   c_dnumda_rad = (cos(dist_rad) - c_r_phi_rad.*cos(s_phi_rad))./(sin(r_phi_rad).*sin(s_phi_rad));
   dnumda_rad = acos(c_dnumda_rad);
   dnumda = rad2deg(dnumda_rad);
   r_numda = s_numda+dnumda;
   if r_numda > 360
       r_numda = r_numda - 360; %��֤������0-360�� 
   end
else
   azi = 360-azi; %����λ�Ǵ���180��ʱ��Ҫ������������λ�øı�
   s_numda_rad = deg2rad(s_numda);
   s_phi_rad = deg2rad(s_phi);
   dist_rad = deg2rad(dist);
   azi_rad = deg2rad(azi);
   c_r_phi_rad = cos(s_phi_rad).*cos(dist_rad) + sin(s_phi_rad).*sin(dist_rad).*cos(azi_rad);
   r_phi_rad = acos(c_r_phi_rad);
   r_phi = rad2deg(r_phi_rad);
   %dnumda�Ǿ��Ȳ�
   c_dnumda_rad = (cos(dist_rad)-c_r_phi_rad.*cos(s_phi_rad))./(sin(r_phi_rad).*sin(s_phi_rad));
   dnumda_rad = acos(c_dnumda_rad);
   dnumda = rad2deg(dnumda_rad);
   r_numda = s_numda-dnumda;
   if r_numda < 0
       r_numda = r_numda + 360; %��֤������0-360��
   end
end
end

