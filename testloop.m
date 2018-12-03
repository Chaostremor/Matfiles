clear; clc;
% vp1 = 1.5: 0.1: 2;
% nvp1 = length(vp1);
% 
% vp3 = 5.9: 0.1: 6.5;
% nvp3 = length(vp3);
% 
% tk1 = 0.1: 0.2: 1.0;
% ntk1 = length(tk1);
% 
% tk2 = 0.5: 0.2: 1.0;
% ntk2 = length(tk2);
% 
% icount = 1;
% aa = zeros(nvp1*nvp3*ntk1, 1);
% for ii = 1:nvp1                                 % thickness of sediment, tk
%     for jj = 1:nvp3                            % source depth, dep
%         for kk = 1:ntk1                       % vp of sediment, vp1
%             for mm = 1:ntk2
%             aa(icount) = tk2(mm)-vp1(ii)+vp3(jj)*tk1(kk);
%             icount = icount+1;
%             end
%         end
%     end
% end
% [val, index] = max(aa(:));
% [itk2, itk1, ivp3, ivp1] = ind2sub([ntk2, ntk1, nvp3, nvp1], index);    % 
% bvp1 = vp1(ivp1);
% bvp3 = vp3(ivp3);
% btk1 = tk1(itk1);
% btk2 = tk2(itk2);
% val1 = -bvp1+bvp3*btk1+btk2;

syms x y
a=2; b=3; c=7; d=1; e=4; f=8;
S1=a*x+b*y-c;
S2=d*x-e*y-f;
S=solve(S1, S2);
x=roundn(double(S.x), -1);
y=roundn(double(S.y), -1);
x
y