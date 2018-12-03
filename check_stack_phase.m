% check_stack_phase
%
% check and verify if the result phase stacked by free_shift_stack is
% reasonable
% eg. Sn, Sg, Sb, Pn, Pb
% 
% Author: C. Song,  2017.4.26
% 

% datadir = 'G:\Alxa\nodecimate\3test400\' ;    % ��������Ŀ¼
% fid1 = fopen(strcat(datadir,'fweight.dat')) ;      % strcat�����ַ�������
% weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊflielist{1}  
% %
% fclose(fid1) ;
% stnm = char(weight{1});     % ̨վ���ǵ�һ��
% [sa, sb] = size(stnm);
% % test to get size of data
% ii=1;
% filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % ��ȡ�ļ�·��
% [ttest, datatest, SAChdrtest] = fget_sac(filename) ;
% disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
% dt = SAChdrtest.times.delta ;
% npts = SAChdrtest.data.trcLen ; 
% %
% [sc, sd] = size(ttest);
% t = zeros(sc,sa);          % t��ʱ�䣬data������
% data = zeros(sc,sa);
% dist = zeros(sa,1);       % dist��ͷ���е����о�
% for ii = 1:sa
%     filename = strcat(datadir,strcat(stnm(ii,:),'.z')) ;  % ��ȡ�ļ�·��
%     [t(:,ii), data(:,ii), SAChdr] = fget_sac(filename) ;   % ��ȡsac�ļ���t��ʱ�䣬data�����ݣ�SAChdr��ͷ�α���
%     dist(ii) = SAChdr.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
% end

figure
for i = 1: sa
%     temp1 = data2(:, i)./(max(abs(data2(:, i)))) .*30.0 ;
%     temp2 = (0.5*(temp1+abs(temp1)));
%     plot(t2(:, i), temp1+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0 0.565 1], 'LineWidth',1); hold on;
%     fill(t2(:, i), temp2+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
% %     plot(t2(figpstart:figpend, i), temp1(figpstart:figpend)+ double(dist(i+index-1)), 'linestyle', '-', 'color', [0.39 0.58 1], 'LineWidth',1); hold on;
% %     fill(t2(figpstart:figpend, i), temp2(figpstart:figpend)+ double(dist(i+index-1)), 'r', 'edgealpha', 0); hold on;
% %fill([t(:, 1), fliplr(t(:, 1))], [y, fliplr(a)], 'r', 'edgealpha', 0);
    plot(t(:,i), data(:,i)./(max(abs(data(:,i)))) .*15.0 + double(dist(i)), 'LineWidth',1);hold on;
end
set(gca, 'XTick', 0:10:200);
set(gca, 'xlim', [0, 200]);
set(gca, 'ylim', [40, 420]);

%% P WAVE PHASE
% distcurve = 390: 0.1: 600;
% p1 = 1/8.2-0.003;         % pPn
% tau1 = 77.35;
% y1 = p1.*(distcurve - 498.6431)+tau1;
% plot(y1, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p2 = 1/8.2-0.004;                   % Pn
% tau2 = 72.15;
% y2 = p2.*(distcurve- 498.6431)+tau2;     
% plot(y2, distcurve, 'linestyle', '--', 'color', 'm', 'LineWidth', 2); hold on;
% 
% % p3 = 1/6.5-0.035;          % ?
% % tau3 = 79.57;
% % y3 = p3.*(distcurve - 498.6431)+tau3;
% % plot( y3, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% p4 = 1/6.5-0.008;          % ?
% tau4 = 81.7;
% y4 = p4.*(distcurve - 498.6431)+tau4;
% plot( y4, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% p5 = 1/6.5-0.034;          % ?
% tau5 = 89.24;
% y5 = p5.*(distcurve - 544.5910)+tau5;
% plot( y5, distcurve, 'linestyle', '--', 'color', 'y', 'LineWidth', 2); hold on;
% 
% a1 = 5.969;
% b1 = 35.05;
% vpg = b1/a1;
% pgcurve = a1.*sqrt((distcurve./b1).^2+1);
% plot(pgcurve, distcurve, 'b-', 'LineWidth', 2); hold on;



% distcurve = 0: 0.1: 600;
% p1 = 1/7.93;         % Pn
% tau1 = 56.87;
% y1 = p1.*(distcurve - 375.5376)+tau1;
% plot(y1, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p2 = 1/7.95;                   % pPn
% tau2 = 62.87;
% y2 = p2.*(distcurve - 375.5376)+tau2;     
% plot(y2, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p3 = 1/6.12;         % ?
% tau3 = 66.64;
% y3 = p3.*(distcurve - 375.5376)+tau3;
% plot( y3, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p4 = 1/6.19;          % ??
% tau4 = 66.73;
% y4 = p4.*(distcurve - 375.5376)+tau4;
% plot( y4, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% p5 = 1/7.18;          % ?
% tau5 = 71.2;
% y5 = p5.*(distcurve - 404.7389)+tau5;
% plot( y5, distcurve, 'linestyle', '--', 'color', 'y', 'LineWidth', 2); hold on;


%% S WAVE PHASE
% distcurve = 200: 0.1: 610;
% p1 = 0.217;         % Sn
% tau1 = 122.8;
% y1 = p1.*(distcurve - 498.6431 )+tau1;
% plot(y1, distcurve, 'linestyle', '--', 'color', 'k', 'LineWidth', 2); hold on;
% 
% p2 = 0.217;                   % sSn
% tau2 = 130.8;
% y2 = p2.*(distcurve- 498.6431)+tau2;     
% plot(y2, distcurve, 'linestyle', '--', 'color', 'y', 'LineWidth', 2); hold on;
% 
% p3 = 0.235;          % Sb
% tau3 = 136.2;
% y3 = p3.*(distcurve - 498.6431)+tau3;
% plot( y3, distcurve, 'linestyle', '--', 'color', 'g', 'LineWidth', 2); hold on;
% 
% a2 = 10.34;              % Sg
% b2 = 35.95;
% vsg = b2/a2;
% sgcurve = a2.*sqrt((distcurve./b2).^2+1);
% plot(sgcurve, distcurve, 'b-', 'LineWidth', 2);
% 
% text(150, 610, 'Sn');
% text(156, 610, 'sSn');
% text(165, 610, 'Sb');
% text(175, 610, 'Sg');


distcurve = 0: 0.1: 400;
p1 = 1/4.53;         % Sn
tau1 = 88.09;
y1 = p1.*(distcurve - 335.1567 )+tau1;
plot(y1, distcurve, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;

p2 = 1/4.13;                   % sSn
tau2 = 95.96;
y2 = p2.*(distcurve - 335.1567)+tau2;     
plot(y2, distcurve, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;
 
p3 = 1/4.52;          % Sb
tau3 = 95.12;
y3 = p3.*(distcurve - 335.1567)+tau3;
plot( y3, distcurve, 'linestyle', '-', 'color', 'k', 'LineWidth', 2); hold on;

a2 = 10.34;              % Sg
b2 = 35.95;
vsg = b2/a2;
sgcurve = a2.*sqrt((distcurve./b2).^2+1);
plot(sgcurve, distcurve, 'b-', 'LineWidth', 2);