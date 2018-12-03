% 进行直角坐标系下和球坐标系下地震矩张量的坐标变换
% Mxx,Myy,Mzz,Mxy,Mxz,Myz是直角坐标系下的矩张量分量
% Mrr,Mtt,Mpp,Mrt,Mrp,Mtp是球坐标系下的矩张量分量
% 这里的x是北方向，y是东方向，z是竖直向下
% 同时，r是径向，t是余纬度指向，p是东经指向，逆时针，从y轴量起
% numda是从y轴逆时针的经度（0-360°），phi是余纬度（0-180°）
%
function [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp] = Coordinate_trans(Mxx,Myy,Mzz,Mxy,Mxz,Myz,numda,phi)
Mcart = [Mxx,Mxy,Mxz;Mxy,Myy,Myz;Mxz,Myz,Mzz]; %Mcart是笛卡尔系地震矩张量
a = deg2rad(numda);
b = deg2rad(phi);
cart2sph = [sin(a)*cos(b) cos(a)*cos(b) sin(b);-sin(a)*sin(b) -cos(a)*sin(b) cos(b);cos(a) -sin(a) 0];
%sph2cart = [sin(a)*cos(b) -sin(a)*sin(b) cos(a);cos(a)*cos(b) -cos(a)*sin(b) -sin(a);sin(b) cos(b) 0];
Msph = cart2sph * Mcart * (cart2sph)'; %Msph是球坐标系地震矩张量
%Mcart = sph2cart * Msph * (sph2cart)';
Mrr = Msph(1,1);
Mtt = Msph(2,2);
Mpp = Msph(3,3);
Mrt = Msph(1,2);
Mrp = Msph(1,3);
Mtp = Msph(2,3);
end