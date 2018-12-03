% test frequency analysis of one signal, need fft
%
% 
% datadir = 'G:\Alxa\3test500\' ;    % 500km��Χδ��������
% fid1 = fopen(strcat(datadir,'fweight.dat')) ;
% weight = textscan(fid1, '%s %f %d %d %d %d %d %f %f \n') ;  % �õ�����һ��cell array,�������Ϊweight{1}  
% %
% fclose(fid1) ;
% stnm = char(weight{1});
% [sa, sb] = size(stnm);
% % test to get size of data
% i=1;
% filename = strcat(datadir,strcat(stnm(i,:),'.t')) ;  % ��ȡ�ļ�·��
% [ttest, datatest, SAChdrtest] = fget_sac(filename) ;
% disttest = SAChdrtest.evsta.dist ;   % ע��ͷ�α����Ĵ洢��ʽ
% dt = SAChdrtest.times.delta ;
% npts = SAChdrtest.data.trcLen ; 
% %
% [sc, sd] = size(ttest);

NFFT = 2^nextpow2(nsnlen+1); % Next power of 2 from length of y
for i=1:nuse
fcoda = fft(dsn(:, i), NFFT);
df = 1/dt/NFFT;
f1 = 0: NFFT/2-1;
f = cat(2, -f1(end: -1: 1), f1).*df; 
amp = abs(fcoda);
%phase = angle(fdata);
amp = fftshift(amp);
%phase = fftshift(phase);
figure (2), plot(f, amp);
set(gca, 'XTICK', -50:2:50);
pause(2);
%figure, plot(f, phase);
end
