% solve equations
%
% Author: C.Song, 2017.4.20

% syms h vp
% a=18;
% x = solve((a-h)/5.87+h/vp==5.969, 2*(a-h)*sqrt((1/5.87)^2-(1/8.33)^2)+2*h*sqrt((1/vp)^2-(1/8.33)^2)==5.25);
% h = vpa(x.h, 4);
% vp = vpa(x.vp, 4);
% h
% vp
a1 = 5.969;
b1 = 35.05;
vpg = b1/a1;
distcurve = 0:0.1:400;
pgcurve = a1.*sqrt((distcurve./b1).^2+1);
figure (1);
plot(pgcurve, distcurve, 'b-'); hold on;

thick = 0.1:0.1:10;
ntk = length(thick);
h = 10:0.1:30;
nh = length(h);
vp = 3.0:0.01:5.5;
nvp = length(vp);
nsample =10000;
pP = linspace(0, 1.0/5.87-0.00005, nsample)';
hflag = 2;
misfit = zeros(ntk, nh, nvp);
for i = 1:ntk                     % thickness of sediment, thick
    for j = 1:nh                  % source depth, h
        for k = 1:nvp           % vp of sediment, vp
            XP=zeros(nsample,1);
            TP=zeros(nsample,1);
            for m = 1:nsample
                 etaP1 = sqrt((1.0/vp(k))^2-pP(m)^2);
                 etaP2 = sqrt((1.0/5.87)^2-pP(m)^2);
                 XP(m) = pP(m) * (thick(i)/etaP1+(h(j)-thick(i))/etaP2);
                 TP(m) = pP(m)*XP(m)+thick(i)*etaP1+(h(j)-thick(i))*etaP2;
            end
            pgint = interp1(XP, TP, distcurve, 'spline');
            misfit(i, j, k) = sum((pgint-pgcurve).^2);
        end
    end
end
[bestfit, index] = min(misfit(:));
[itk, jh, kvp] = ind2sub([ntk, nh, nvp], index);
besttk = thick(itk);
besth = h(jh);
bestvp = vp(kvp);
XP=zeros(nsample,1);
TP=zeros(nsample,1);
for m = 1:nsample
    etaP1 = sqrt((1.0/bestvp)^2-pP(m)^2);
    etaP2 = sqrt((1.0/5.87)^2-pP(m)^2);
    XP(m) = pP(m) * (besttk/etaP1+(besth-besttk)/etaP2);
    TP(m) = pP(m)*XP(m)+besttk*etaP1+(besth-besttk)*etaP2;
end
pgint = interp1(XP, TP, distcurve, 'spline');
figure(2)
plot(pgcurve, distcurve, 'b-'); hold on;
plot(pgint, distcurve, 'r-'); 