dt=0.001;
fm=20;%ÕðÔ´Ö÷Æµ
t=-0.06:dt:0.06;
nn=length(t);
A=1;
w=A*(1-2*(pi*fm*t).^2).*exp(-(pi*fm*t).^2);
wave=w;
ft=wave;
FW=fft(ft);
figure (1);
subplot(2,2,1);
plot(t,ft),title('ft');
subplot(2,2,2);
plot(FW),title('FW');
amp=abs(FW);
nn=1:121;
%%
dw=(nn-61)/(0.001*120);
%%
subplot(2,2,3);
%%
plot(dw,fftshift(amp)),title('amplitude');axis([-50,50,0,20])
%%
phrase=angle(FW)*180/pi;
subplot(2,2,4);
plot(dw,phrase),title('phrase');
% ÂË²¨Æ÷
order = 3 ;
delta = 0.001 ;
lowf = 5;
highf = 25 ;
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
%%
plot(dw,fftshift(eamp)),title('eamplitude');axis([-50,50,0,20])
%%
ephrase=angle(eFW)*180/pi;
subplot(2,2,4);
plot(dw,ephrase),title('ephrase');