%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 将某深度震源在若干个震中距上计算的不含方位角信息的地震图按方位角投影恢复 %
% Revised history:                                                                                             %
%     2016.11.16         initial coding , Gezx                                                          %
%     2016.11.17         modification , Sc                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
clear
% 震源位置eq location
lateq = 0.0 ;
loneq = 0.0 ;

% 台站位置st location
latst = [60.0; 48.5904; 25.6589; 0.0] ;
lonst = [0.0; 40.8934; 56.3099; 60.0] ;

% 计算每个台站的震中距dis，方位角azi，自带函数distance
[dis1, azi] = distance(lateq, loneq, latst, lonst) ;
% qpwvinf.f 定义的方位角azi从南向东逆时针计算
% 与传统的azi从北顺时针不同，两者相加等于180度
azi = (180-azi)*pi/180 ;   
ssa = sin(azi) ;
csa = cos(azi) ;
ss2a = sin(2*azi) ;
cs2a = cos(2*azi) ;

% 计算每个台站的反方位角bazi（从台站到震源）
[dis2, bazi] = distance(latst, lonst, lateq, loneq) ;
% qpwvinf.f 定义的反方位角bazi也和传统的不一致
% 代码中的是传统定义减去180度
bazi = (bazi-180)*pi/180 ;
ssb = sin(bazi) ;
csb = cos(bazi) ;

% read in green's functions for 6 basic moment tensors
% to avoid misunderstanding , recommend qssp ouput be TPR components 

% 1: expl , corresponding tensor=[1 0 0;0 1 0;0 0 1] in TPR coordinate
expl = importdata('G:\qssp2010-code+input\Decomposition\expl\expl.up') ;   
explup = expl.data(:, 2:end) ;
expl = importdata('G:\qssp2010-code+input\Decomposition\expl\expl.ut') ;
explut = expl.data(:, 2:end) ;
expl = importdata('G:\qssp2010-code+input\Decomposition\expl\expl.ur') ;
explur = expl.data(:, 2:end) ;

% 2: clvd , corresponding tensor=[-0.5 0 0;0 -0.5 0;0 0 1]
clvd = importdata('G:\qssp2010-code+input\Decomposition\clvd\clvd.up') ;
clvdup = clvd.data(:, 2:end) ;
clvd = importdata('G:\qssp2010-code+input\Decomposition\clvd\clvd.ut') ;
clvdut = clvd.data(:, 2:end) ;
clvd = importdata('G:\qssp2010-code+input\Decomposition\clvd\clvd.ur') ;
clvdur = clvd.data(:, 2:end) ;

% 3: s12 , corresponding tensor=[0 1 0;1 0 0;0 0 0]
ss12 = importdata('G:\qssp2010-code+input\Decomposition\ss12\ss12.up') ;
ss12up = ss12.data(:, 2:end) ;
ss12 = importdata('G:\qssp2010-code+input\Decomposition\ss12\ss12.ut') ;
ss12ut = ss12.data(:, 2:end) ;
ss12 = importdata('G:\qssp2010-code+input\Decomposition\ss12\ss12.ur') ;
ss12ur = ss12.data(:, 2:end) ;

% 4: s11 , corresponding tensor=[1 0 0;0 -1 0;0 0 0]
ss11 = importdata('G:\qssp2010-code+input\Decomposition\ss11\ss11.up') ;
ss11up = ss11.data(:, 2:end) ;
ss11 = importdata('G:\qssp2010-code+input\Decomposition\ss11\ss11.ut') ;
ss11ut = ss11.data(:, 2:end) ;
ss11 = importdata('G:\qssp2010-code+input\Decomposition\ss11\ss11.ur') ;
ss11ur = ss11.data(:, 2:end) ;

% 5: ds31，corresponding tensor=[0 0 1;0 0 0;1 0 0]
ds31 = importdata('G:\qssp2010-code+input\Decomposition\ds31\ds31.up') ;
ds31up = ds31.data(:, 2:end) ;
ds31 = importdata('G:\qssp2010-code+input\Decomposition\ds31\ds31.ut') ;
ds31ut = ds31.data(:, 2:end) ;
ds31 = importdata('G:\qssp2010-code+input\Decomposition\ds31\ds31.ur') ;
ds31ur = ds31.data(:, 2:end) ;

% 6: ds23 , corresponding tensor=[0 0 0;0 0 1;0 1 0] 
ds23 = importdata('G:\qssp2010-code+input\Decomposition\ds23\ds23.up') ;
ds23up = ds23.data(:, 2:end) ;
ds23 = importdata('G:\qssp2010-code+input\Decomposition\ds23\ds23.ut') ;
ds23ut = ds23.data(:, 2:end) ;
ds23 = importdata('G:\qssp2010-code+input\Decomposition\ds23\ds23.ur') ;
ds23ur = ds23.data(:, 2:end) ;

% to compare the waveform, compute whole tensor as well
allkind = importdata('G:\qssp2010-code+input\Decomposition\all\all.uz') ;
alluz = allkind.data(:, 2:end) ;
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moment tenser decomposition                                    %
% M11~M33 actually mean 1-T, 2-P, 3-R                         %
% all kind of mt can be divided into these 6 basic ones    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M0, M11~M33 can vary in real cases
M0 = 1.0 ;
M11 = 1.0 ;
M22 = 1.0 ;
M33 = 1.0 ;
M12 = 1.0 ;
M13 = 1.0 ;
M23 = 1.0 ;
cexpl = (M11+M22+M33)/3 ;
cclvd = M33-cexpl ;
css11 = (M11-M22)/2 ;
css12 = M12 ;
cds31 = M13 ;
cds23 = M23 ;

% expl and clvd are independent on azimuth 
explup = explup*cexpl ;
explut = explut*cexpl ;
explur = explur*cexpl ;
clvdup = clvdup*cclvd ;
clvdut = clvdut*cclvd ;
clvdur = clvdur*cclvd ;

% choose any data (same size) to get the size
[nrec, nst] = size(ds23ur) ;

% initial
ss12upi = zeros(nrec,nst) ;
ss12uri = zeros(nrec,nst) ;
ss12uti = zeros(nrec,nst) ;
ds31upi = zeros(nrec,nst) ;
ds31uri = zeros(nrec,nst) ;
ds31uti = zeros(nrec,nst) ;
ss11uri = zeros(nrec,nst) ;
ss11upi = zeros(nrec,nst) ;
ss11uti = zeros(nrec,nst) ;
ds23uri = zeros(nrec,nst) ;
ds23upi = zeros(nrec,nst) ;
ds23uti = zeros(nrec,nst) ;
ur = 0 ; 
up = 0 ;
ut = 0 ; 
uz = zeros(nrec,nst) ;
ux = zeros(nrec,nst) ;
uy = zeros(nrec,nst) ;
% loop
for i=1:4
    % ss12
    ss12upi(:, i) = ss12up(:, i).*cs2a(i)*css12 ;
    ss12uri(:, i) = ss12ur(:, i).*ss2a(i)*css12 ;
    ss12uti(:, i) = ss12ut(:, i).*ss2a(i)*css12 ;
    
    % ds31
    ds31upi(:, i) = ds31up(:, i).*ssa(i)*cds31 ;
    ds31uri(:, i) = ds31ur(:, i).*csa(i)*cds31 ;
    ds31uti(:, i) = ds31ut(:, i).*csa(i)*cds31 ;
    
    % ss11
    ss11uri(:, i) = ss11ur(:, i).*cs2a(i)*css11 ;
    ss11upi(:, i) = -ss11up(:, i).*ss2a(i)*css11 ;
    ss11uti(:, i) = ss11ut(:, i).*cs2a(i)*css11 ;
    
    % ds23
    ds23uri(:, i) = ds23ur(:, i).*ssa(i)*cds23 ;
    ds23upi(:, i) = ds23up(:, i).*ssa(i)*cds23 ;
    ds23uti(:, i) = -ds23ut(:, i).*csa(i)*cds23 ;
    
    % compound
    ur = explur(:, i)+clvdur(:, i)+ss12uri(:, i)+ds31uri(:, i)+ss11uri(:, i)+ds23uri(:, i) ;
    up = explup(:, i)+clvdup(:, i)+ss12upi(:, i)+ds31upi(:, i)+ss11upi(:, i)+ds23upi(:, i) ;
    ut = explut(:, i)+clvdut(:, i)+ss12uti(:, i)+ds31uti(:, i)+ss11uti(:, i)+ds23uti(:, i) ;
    
    % project to coordinate NED
    uz(:, i) = ur ;
    ux(:, i) = ut.*csb(i)+up.*ssb(i) ;
    uy(:, i) = ut.*ssb(i)-up.*csb(i) ;
    
    % plot to find whether the result is ture 
    figure
    plot(allkind.data(:,1), alluz(:, i), 'g.') ;
    hold on;
    plot(allkind.data(:,1), uz(:,i), 'r') ;
end



