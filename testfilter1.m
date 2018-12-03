t=pi/2:0.01:pi;
ft=sin(t);
FW=fft(ft);
figure (1);
subplot(2,2,1);
plot(t,ft),title('ft');
subplot(2,2,2);
plot(FW),title('FW');
amp=abs(FW);
subplot(2,2,3);
plot(amp),title('amplitude');axis([0,5,0,100])
phrase=angle(FW)*180/pi;
subplot(2,2,4);
plot(phrase),title('phrase');
% ÂË²¨Æ÷
order = 3 ;
delta = 0.01 ;
lowf = 0.1;
highf = 0.2 ;
fs = 1.0/delta ;
nyq = fs/2.0 ;
[B,A] = butter( order, [lowf/nyq, highf/nyq] ) ;
eft=filtfilt(B,A,ft);
%
eFW=fft(eft);
figure (2);
subplot(2,2,1);
plot(t,eft),title('eft');
subplot(2,2,2);
plot(eFW),title('eFW');
eamp=abs(eFW);
subplot(2,2,3);
plot(eamp),title('eamplitude');%axis([0,5,0,20])
ephrase=angle(eFW)*180/pi;
subplot(2,2,4);
plot(ephrase),title('ephrase');