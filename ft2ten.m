% ���ϲ�Ĳ���ת��Ϊ�������ĸ�����
% ����ϲ����������phis�����delta��������lamda
% �����������Mxx��Myy��Mzz��Mxz��Mxy��Myz
% ���������Լ����x�Ǳ�����y�Ƕ�����z����ֱ����
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