% 简正振型中的地震矩张量表达
% 从笛卡尔坐标系转换到球坐标系（都是对于震源的坐标系，指向之中在变化）
% 笛卡尔坐标系的定义：x的方向是北；y的方向是东；z的方向是向下
% 球坐标系的定义：r是径向；t是余纬度方向；p是东经度方向，逆时针
%
function [Msph] = Normal_mode_mtensor(Mcart)
Mxx = Mcart(1,1); %Mcart是笛卡尔系地震矩张量
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
Msph = [Mrr Mrt Mrp;Mrt Mtt Mtp;Mrp Mtp Mpp]; %Msph是球坐标系地震矩张量
end