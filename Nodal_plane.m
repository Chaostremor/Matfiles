% relations between the fault plane, the auxiliary plane, and the stress
% axes
%
% AUTHOR: 
%     C Song, 2017.7.7
% REFERENCE: 
%     An introduction to Seismology, Earthquakes, and Earth Structure

%
function [stk2, dip2, rake2]=Nodal_plane(stk1, dip1, rake1)
% clear;
% stk1 = 356;
% dip1 = 89;
% rake1 = -167;
% convert deg to rad
stk1 = deg2rad(stk1);
dip1 = deg2rad(dip1);
rake1 = deg2rad(rake1);
% fault normal vector
n1(1,1) = -sin(dip1)*sin(stk1);
n1(2,1) = -sin(dip1)*cos(stk1); 
n1(3,1) = cos(dip1);
% fault slip vector
d1(1,1) = cos(rake1)*cos(stk1)+sin(rake1)*cos(dip1)*sin(stk1);
d1(2,1) = -cos(rake1)*sin(stk1)+sin(rake1)*cos(dip1)*cos(stk1);
d1(3,1) = sin(rake1)*sin(dip1);
% null axis vector
b(1,1) = -sin(rake1)*cos(stk1)+cos(rake1)*cos(dip1)*sin(stk1);
b(2,1) = sin(rake1)*sin(stk1)+cos(rake1)*cos(dip1)*cos(stk1);
b(3,1) = cos(rake1)*sin(dip1);
% T axis, minimum compressive stress in compressional quadrants
t = n1+d1;
% P axis, maximum compressive stress in dilatational quadrants
p = n1-d1;

% get the other nodal plane from one
% first get dip
% dip2 = acos(sin(rake1)*sin(dip1));
% dip2 = rad2deg(dip2);
% if dip2>90 && dip2<180
%     dip2 = 180-dip2;
% end
% 
% dip2 = deg2rad(dip2); 
% srake2 = cos(dip1)/sin(dip2);
% crake2 = -sin(dip1)*cos(rake1)/sin(dip2);
% rake2 = acos(crake2);
% rake2 = rad2deg(rake2);
% if srake2 <0 
%     rake2 = 360-rake2;
% end
% if rake2 > 180
%     rake2 = rake2-360;
% end
% 
% sstk12 = cos(dip1)/sin(dip2);
% cstk12 = -1/tan(dip1)/tan(dip2);
% stk12 = acos(cstk12);
% stk12 = rad2deg(stk12);
% if sstk12 < 0 
%     stk12 = 360-stk12;
% end
% if rad2deg(stk1) >=180   
%     stk2 = rad2deg(stk1)-stk12;
% else
%     stk2 = rad2deg(stk1)-stk12 +180;
% end


dip2 = acos(sin(rake1)*sin(dip1)); 

srake2 = cos(dip1)/sin(dip2);    % >0
crake2 = -sin(dip1)*cos(rake1)/sin(dip2);
rake2 = acos(crake2);
rake2 = rad2deg(rake2);

sstk12 = cos(rake1)/sin(dip2);
cstk12 = -1/tan(dip1)/tan(dip2);
stk12 = acos(cstk12);
stk12 = rad2deg(stk12);
if sstk12 < 0 
    stk12 = 360-stk12;
end
stk2 =rad2deg(stk1) - stk12;

dip2 = rad2deg(dip2);

if dip2>=90 && dip2<180
    dip2 = 180-dip2;
    rake2 = 360-rake2;
    stk2 = 180+stk2;
end

if rake2 > 180                    % to make rake [-180, 180]
    rake2 = rake2-360;
end

if stk2 < 0
    stk2 = stk2+360;
end

end