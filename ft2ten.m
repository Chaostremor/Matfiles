% 将断层的参数转化为矩张量的各分量
% 输入断层参数：走向phis，倾角delta，滑动角lamda
% 输出矩张量：Mxx，Myy，Mzz，Mxz，Mxy，Myz
% 这里的坐标约定是x是北方向，y是东方向，z是竖直向下
%
function [M] = ft2ten(phis, delta, lamda, M0)
 x = deg2rad(delta) ;
 y = deg2rad(lamda) ; 
 z = deg2rad(phis) ;
 sinx = sin(x) ;
 cosx = cos(x) ;
 sin2x = sin(2*x) ; 
 cos2x = cos(2*x) ;
 siny = sin(y) ;
 cosy = cos(y) ;
 sinz = sin(z) ;
 cosz = cos(z) ;
 sin2z = sin(2*z) ;
 cos2z = cos(2*z) ;
 sinz2 = sinz^2.0 ;
 cosz2 = cosz^2.0 ;
 
 Mxx = -M0*(sinx*cosy*sin2z+sin2x*siny*sinz2) ; 
 Mxy = +M0*(sinx*cosy*cos2z+0.5*sin2x*siny*sin2z) ; 
 Mxz = -M0*(cosx*cosy*cosz+cos2x*siny*sinz) ;
 Myy = +M0*(sinx*cosy*sin2z-sin2x*siny*cosz2) ;
 Myz = -M0*(cosx*cosy*sinz-cos2x*siny*cosz) ;
 Mzz = +M0*sin2x*siny ;
 
 M = [Mxx Mxy Mxz;Mxy Myy Myz;Mxz Myz Mzz];
 
end