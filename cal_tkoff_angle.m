% 计算离源角ih
%
% r0是地球半径
% h是震源所在的深度，rh=r0-h
% v0是地面处的地震波波速（根据波的种类）
% vh是震源处的地震波深度
% i0是入射角
% 球对称介质的Snell定律是：rh*sin(ih)/vh=r0*sin(i0)/v0
% 入射角可由地震波的时距曲线t（delta，h）求得：
% sin(i0)=v0*dt/ddelta，其中delta以长度例如km为单位
% 则最终可以得到：
% sin(ih)=r0/(r0-h)*vh*dt/ddelta
function [ih] = cal_tkoff_angle(r0,h,vh,dtddelta)
sinih = r0/(r0-h)*vh*dtddelta;
ih = asin(sinih);
ih = rad2deg(ih);
return
end