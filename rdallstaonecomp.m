% ��������̨վ��ĳһ��λ�Ʒ���
%
clear; close all;
%
%k = 1;
%
% un = importdata('G:\qssp2010-code+input\Japantest\3Japan_out.un');
% Tun = un.data(:,1);
% %
% Tun = Tun - 148 + 460;
% %
% wavedataun = un.data(:,2:end);
% [mn,nn] = size(wavedataun);
% unst1 = wavedataun(:,k);
% %
% ue = importdata('G:\qssp2010-code+input\Japantest\3Japan_out.ue');
% Tue = ue.data(:,1);
% %
% Tue = Tue - 148 + 460;
% %
% wavedataue = ue.data(:,2:end);
% [me,ne] = size(wavedataue);
% uest1 = wavedataue(:,k);
% %
uz = importdata('G:\qssp2010-code+input\Japantest\4Japan_out.uz');
stainfo = load('G:\qssp2010-code+input\Japantest\38stainfo.txt');
Tuz = uz.data(:,1);
wavedatauz = uz.data(:,2:end);
[mz,nz] = size(wavedatauz);
treduction = stainfo(:,1);
dist = stainfo(:,2);
% for k = 1:nz
%    uzk = wavedatauz(:,k);
%    tuzk = Tuz -148 + treduction(k);
%    subplot(nz,1,k);
%    plot(tuzk,uzk,'-b');
%    set(gca,'XTick',0:50:2000);
%    set(gca,'XLim',[0 2000]);
% end
for k = 1:nz
    uzk = wavedatauz(:,k)./(max(abs(wavedatauz(:,k)))) + dist(k)/2;
    tuzk = Tuz -148 + treduction(k);
    % 148�ǲ���������Դ��T_origin����ʽ�����ʱ�����Ը�����Դ��ʼʱ��Ϊ�ο�
    % �����ʱ���׼����Ҫע�⡣
    plot(tuzk,uzk,'-b');hold on;
end
set(gca,'XTick',0:50:2000);
set(gca,'XLim',[0 2000]);