% 绘制Azimuth中的所有矩张量分量产生的理论位移的同一个位移分量，比如uz
% 并对数据进行带通滤波
%
clear; close all;
%
%
mppuz = importdata('G:\qssp2010-code+input\Azimuth9\mpp\azi_out_mpp.uz');
mrpuz = importdata('G:\qssp2010-code+input\Azimuth9\mrp\azi_out_mrp.uz');
mrruz = importdata('G:\qssp2010-code+input\Azimuth9\mrr\azi_out_mrr.uz');
mrtuz = importdata('G:\qssp2010-code+input\Azimuth9\mrt\azi_out_mrt.uz');
mtpuz = importdata('G:\qssp2010-code+input\Azimuth9\mtp\azi_out_mtp.uz');
mttuz = importdata('G:\qssp2010-code+input\Azimuth9\mtt\azi_out_mtt.uz');
% muz = importdata('G:\qssp2010-code+input\Azimuth8\m\azi_out_m.uz');

% mppuz = importdata('G:\qssp2010-code+input\Azimuth1\mpp\azi_out_mpp.uz');
% mrpuz = importdata('G:\qssp2010-code+input\Azimuth1\mrp\azi_out_mrp.uz');
% mrruz = importdata('G:\qssp2010-code+input\Azimuth1\mrr\azi_out_mrr.uz');
% mrtuz = importdata('G:\qssp2010-code+input\Azimuth1\mrt\azi_out_mrt.uz');
% mtpuz = importdata('G:\qssp2010-code+input\Azimuth1\mtp\azi_out_mtp.uz');
% mttuz = importdata('G:\qssp2010-code+input\Azimuth1\mtt\azi_out_mtt.uz');

% mppuz = importdata('G:\qssp2010-code+input\Azimuth2\mpp\azi_out_mpp.uz');
% mrpuz = importdata('G:\qssp2010-code+input\Azimuth2\mrp\azi_out_mrp.uz');
% mrruz = importdata('G:\qssp2010-code+input\Azimuth2\mrr\azi_out_mrr.uz');
% mrtuz = importdata('G:\qssp2010-code+input\Azimuth2\mrt\azi_out_mrt.uz');
% mtpuz = importdata('G:\qssp2010-code+input\Azimuth2\mtp\azi_out_mtp.uz');
% mttuz = importdata('G:\qssp2010-code+input\Azimuth2\mtt\azi_out_mtt.uz');
Tun = mppuz.data(:,1);
%
%
wavedatauz(1,:,:) = mppuz.data(:,2:end);
wavedatauz(2,:,:) = mrpuz.data(:,2:end);
wavedatauz(3,:,:) = mrruz.data(:,2:end);
wavedatauz(4,:,:) = mrtuz.data(:,2:end);
wavedatauz(5,:,:) = mtpuz.data(:,2:end);
wavedatauz(6,:,:) = mttuz.data(:,2:end);
% wavedatauz(7,:,:) = muz.data(:,2:end);
wavedatauz(8,:,:) = wavedatauz(1,:,:)+wavedatauz(3,:,:)+wavedatauz(6,:,:);
[ms,ns,ps] = size(wavedatauz);
%
% FW = fft(wavedatauz(2,1:300,1));
% amp = abs(FW);
% figure (1);
% plot(amp),title('amplitude');
%
order = 3 ;
delta = 0.1 ;
lowf = 0.5;
highf = 2.0 ;
fs = 1.0/delta ;
nyq = fs/2.0 ;
[B,A] = butter( order, [lowf/nyq, highf/nyq] ) ;
for i=1:ms
    for j=1:ps
        wavedatauz(i,:,j) = filtfilt(B , A, wavedatauz(i,:,j)) ;
    end
end
%
figure (2);
for i = 1:ms
    tuz = Tun + 554;
    % 冒号表达式构造时间向量，从Treduction开始，长度为Twin+T_origin_max
    uzplot = wavedatauz(i,:,1)./(max(abs(wavedatauz(i,:,1))))+ i;
    plot(tuz,uzplot);hold on;
end
set(gca,'XTick',0:50:1200);
set(gca,'XLim',[500 1100]);  % 设置x轴的间隔，长度
% title('慢度积分格林函数0.1，输出0.1，未滤波');
% title('慢度积分格林函数0.18，输出0.18，0.01-0.20滤波');
% title('慢度积分格林函数0.18，输出0.1，0.01-0.20滤波');

%tt=tauptime('mod','iasp91','ph','P','evt',[0 0],'sta',[60 0],'h',14);
