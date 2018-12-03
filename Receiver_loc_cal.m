%已知一个地震震中的经纬度，以及台站到震中的方位角和震中距
%要求台站的经纬度
%利用球面三角的余弦公式计算
%设震中的经纬度分别为（s_numda，s_phi），（0-360，0-180）
%震中距为dist(0-180)，
%方位角为azi(0-360)，
%待求的台站的经纬度分别为（r_numda,r_phi），（0-360，0-180）
%如果采用的各角度的取值范围不同，需要做简单的修改
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
   %dnumda是经度差,算法限制，s_phi不能是0
   c_dnumda_rad = (cos(dist_rad) - c_r_phi_rad.*cos(s_phi_rad))./(sin(r_phi_rad).*sin(s_phi_rad));
   dnumda_rad = acos(c_dnumda_rad);
   dnumda = rad2deg(dnumda_rad);
   r_numda = s_numda+dnumda;
   if r_numda > 360
       r_numda = r_numda - 360; %保证经度是0-360度 
   end
else
   azi = 360-azi; %当方位角大于180度时，要求解的球面三角位置改变
   s_numda_rad = deg2rad(s_numda);
   s_phi_rad = deg2rad(s_phi);
   dist_rad = deg2rad(dist);
   azi_rad = deg2rad(azi);
   c_r_phi_rad = cos(s_phi_rad).*cos(dist_rad) + sin(s_phi_rad).*sin(dist_rad).*cos(azi_rad);
   r_phi_rad = acos(c_r_phi_rad);
   r_phi = rad2deg(r_phi_rad);
   %dnumda是经度差
   c_dnumda_rad = (cos(dist_rad)-c_r_phi_rad.*cos(s_phi_rad))./(sin(r_phi_rad).*sin(s_phi_rad));
   dnumda_rad = acos(c_dnumda_rad);
   dnumda = rad2deg(dnumda_rad);
   r_numda = s_numda-dnumda;
   if r_numda < 0
       r_numda = r_numda + 360; %保证经度是0-360度
   end
end
end

