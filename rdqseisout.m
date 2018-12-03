clear;close all;
% 读入理论地震图计算程序的输出文件
%A = importdata('G:\qseis-06-examples\compare_qseis\qseisout.tz');
%TA = A.data(:,1);
%wavedataA = A.data(:,2:end);
B = importdata('G:\qssp2010-code+input\Japan\Japan_out.uz');
TB = B.data(:,1);
wavedataB = B.data(:,2:end);
[m,n] = size(wavedataB);
%C = importdata('G:\qssp2010-code+input\compare_qssp\qsspout1.uz');
%TC = C.data(:,1);
%wavedataC = (-1).*C.data(:,2:end);
figure (1)
%normwave(:,1) = wavedata(:,1)./(max(wavedata(:,1))) + 1-1;
% for i = 1:10
%     % 用振幅的最大值对整个波形做归一化
%     normwaveA(:,i) = wavedataA(:,i)./(max(wavedataA(:,i))) + i; 
%     plot(TA,normwaveA(:,i),'b-');hold on;
% end
for i= 1:n
   normwaveB(:,i) = wavedataB(:,i)./(max(abs(wavedataB(:,i)))) + i;
   plot(TB,normwaveB(:,i),'b-');hold on;
end
% for i= 1:10
%     normwaveC(:,i) = wavedataC(:,i)./(max(wavedataC(:,i))) + i;
%     plot(TC,normwaveC(:,i),'r-');hold on;
% end
xlabel('time   /s');
%axis([0,1000,0,11]);
ylabel('station number');

