function [T, X]=cal_dist_time(n, k, h, thick, v, type)
% calculate the distance and travel time of direct P and Pn
% author: C. Song,   2017.3.5
%
% model like below :
%--------------------------------------
%--------------------------------------
%     **
%--------------------------------------
%--------------------------------------
%
% NOTES
%      n ------> number of total layers (included in calculation)
%      k ------> number of layer of source
%      h ------> depth of source
%      thick ----> thickness of each layer
%      v ------> velocity of each layer
%      type -----> phase type, 'direct', 'head', 'ref_head'
%      X ------> distance
%      T ------>  travel time


% clear;close all;
% n=4;
% vs = [3.04; 3.51; 3.77; 4.756];
% vp = [5.09208; 5.90612; 6.35644; 8.06419];
% thick = [5.1; 20.0; 15.0; 0.0];
% k = 2;
% h = 16;

nsample =10000;

if (strcmp(type, 'direct') ==1)
    % for Pg/Sg
    ray = linspace(0, 1.0/v(k)-0.0000001, nsample)';
    if k==1 
        thicksum = 0;
    else
        thicksum = 0;
        for i = 1:k-1
            thicksum = thicksum+thick(i);
        end
    end
    X=zeros(nsample,1);
    T=zeros(nsample,1);
    for i = 1:nsample
        X(i) = 0;
        hi=zeros(k,1);
        eta=zeros(k,1);
        for j = 1:k
            hi(j)=thick(j);
            eta(j) = sqrt((1.0/v(j))^2-ray(i)^2);
            if (j == k)
                hi(j) = h-thicksum;
            end
            X(i) = X(i)+ray(i)*hi(j)/eta(j);
        end
        T(i) = ray(i)*X(i);
        for j = 1:k
            T(i) = T(i)+hi(j)*eta(j);
        end
    end
    
elseif (strcmp(type, 'head') ==1)
        
        % for Pn/Sn, Pb/Sb
        ray = 1/v(n);
        thicksum = 0;
        for i = 1:k
            thicksum = thicksum+thick(i);
        end
        dX = linspace(0, 400, nsample)';
        Xmin = 0;
        hi=zeros(n-1,1);
        eta=zeros(n-1,1);
        for i = 1:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            Xmin = Xmin+ray*hi(i)/eta(i);
        end
        hi=zeros(n-k,1);
        eta=zeros(n-k,1);
        for i = k:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            if (i == k)
                hi(i) = thicksum-h;
            end
            Xmin = Xmin+ray*hi(i)/eta(i);
        end
        X = Xmin + dX;

        T = ray*Xmin;
        hi=zeros(n-1,1);
        eta=zeros(n-1,1);
        for i = 1:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            T = T+hi(i)*eta(i);
        end
        hi=zeros(n-k,1);
        eta=zeros(n-k,1);
        for i = k:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            if (i == k)
                hi(i) = thicksum-h;
            end
            T = T+hi(i)*eta(i);
        end
        T = T + ray*dX;
        
elseif (strcmp(type, 'ref_head') ==1)
        
        % for sSn, reflected by ground surface, then like Sn 
        ray = 1/v(n);
        if k==1 
            thicksum = 0;
        else
            thicksum = 0;
            for i = 1:k-1
                thicksum = thicksum+thick(i);
            end
        end
        dX = linspace(0, 400, nsample)';
        Xmin = 0;
        hi=zeros(n-1,1);
        eta=zeros(n-1,1);
        for i = 1:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            Xmin = Xmin+2*ray*hi(i)/eta(i);
        end
        hi=zeros(k,1);
        eta=zeros(k,1);
        for i = 1:k
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            if (i == k)
                hi(i) = h-thicksum;
            end            
            Xmin = Xmin+ray*hi(i)/eta(i);
        end
        X = Xmin + dX;

        T = ray*Xmin;
        hi=zeros(n-1,1);
        eta=zeros(n-1,1);
        for i = 1:n-1
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            T = T+2*hi(i)*eta(i);
        end
        hi=zeros(k,1);
        eta=zeros(k,1);
        for i = 1:k
            hi(i)=thick(i);
            eta(i) = sqrt((1.0/v(i))^2-ray^2);
            if (i == k)
                hi(i) = h-thicksum;
            end
            T = T+hi(i)*eta(i);
        end
        T = T + ray*dX;
        
end



% % for direct S, Sg
% pS = linspace(0, 1.0/vs(k)-0.0000005, nsample)';
% if k==1 
%     thicksumS = 0;
% else
%     thicksumS = 0;
%     for i = 1:k-1
%         thicksumS = thicksumS+thick(i);
%     end
% end
% XS=zeros(nsample,1);
% TS=zeros(nsample,1);
% for i = 1:nsample
%     XS(i) = 0;
%     hiS=zeros(k,1);
%     etaS=zeros(k,1);
%     for j = 1:k
%         hiS(j)=thick(j);
%         etaS(j) = sqrt((1.0/vs(j))^2-pS(i)^2);
%         if (j == k)
%             hiS(j) = h-thicksumS;
%         end
%         XS(i) = XS(i)+pS(i)*hiS(j)/etaS(j);
%     end
%     TS(i) = pS(i)*XS(i);
%     for j = 1:k
%         TS(i) = TS(i)+hiS(j)*etaS(j);
%     end
% end
% 
% % for Sn wave
% pSn = 1/vs(n);
% thicksumSn = 0;
% for i = 1:k
%     thicksumSn = thicksumSn+thick(i);
% end
% dXSn = linspace(0, 400, nsample)';
% XSnmin = 0;
% hiSn=zeros(n-1,1);
% etaSn=zeros(n-1,1);
% for i = 1:n-1
%      hiSn(i)=thick(i);
%      etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
%      XSnmin = XSnmin+pSn*hiSn(i)/etaSn(i);
% end
% hiSn=zeros(n-k,1);
% etaSn=zeros(n-k,1);
% for i = k:n-1
%      hiSn(i)=thick(i);
%      etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
%      if (i == k)
%          hiSn(i) = thicksumSn-h;
%      end
%      XSnmin = XSnmin+pSn*hiSn(i)/etaSn(i);
% end
% XSn = XSnmin + dXSn;
% 
% TSn = pSn*XSnmin;
% hiSn=zeros(n-1,1);
% etaSn=zeros(n-1,1);
% for i = 1:n-1
%      hiSn(i)=thick(i);
%      etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
%      TSn = TSn+hiSn(i)*etaSn(i);
% end
% hiSn=zeros(n-k,1);
% etaSn=zeros(n-k,1);
% for i = k:n-1
%      hiSn(i)=thick(i);
%      etaSn(i) = sqrt((1.0/vs(i))^2-pSn^2);
%      if (i == k)
%          hiSn(i) = thicksumSn-h;
%      end
%       TSn = TSn+hiSn(i)*etaSn(i);
% end
% TSn = TSn + pSn*dXSn;

% % head wave of interface between  upper crust and lower crust, Conrad, Pb, Sb
% pSb = 1/vs(n-1);
% thicksumSb = 0;
% for i = 1:k
%     thicksumSb = thicksumSb+thick(i);
% end
% dXSb = linspace(0, 400, nsample)';
% XSbmin = 0;
% hiSb=zeros(n-2,1);
% etaSb=zeros(n-2,1);
% for i = 1:n-2
%      hiSb(i)=thick(i);
%      etaSb(i) = sqrt((1.0/vs(i))^2-pSb^2);
%      XSbmin = XSbmin+pSb*hiSb(i)/etaSb(i);
% end
% hiSb=zeros(n-k-1,1);
% etaSb=zeros(n-k-1,1);
% for i = k:n-2
%      hiSb(i)=thick(i);
%      etaSb(i) = sqrt((1.0/vs(i))^2-pSb^2);
%      if (i == k)
%          hiSb(i) = thicksumSb-h;
%      end
%      XSbmin = XSbmin+pSb*hiSb(i)/etaSb(i);
% end
% XSb = XSbmin + dXSb;
% 
% TSb = pSb*XSbmin;
% hiSb=zeros(n-2,1);
% etaSb=zeros(n-2,1);
% for i = 1:n-2
%      hiSb(i)=thick(i);
%      etaSb(i) = sqrt((1.0/vs(i))^2-pSb^2);
%      TSb = TSb+hiSb(i)*etaSb(i);
% end
% hiSb=zeros(n-k-1,1);
% etaSb=zeros(n-k-1,1);
% for i = k:n-2
%      hiSb(i)=thick(i);
%      etaSb(i) = sqrt((1.0/vs(i))^2-pSb^2);
%      if (i == k)
%          hiSb(i) = thicksumSb-h;
%      end
%       TSb = TSb+hiSb(i)*etaSb(i);
% end
% TSb = TSb + pSb*dXSb;



% plot(TP, XP);
% fid = fopen('C:\Users\Song Chao\Desktop\fweight_select.dat');
% weight = textscan(fid, '%s %f %d %d %d %d %d %f %f');
% Xinpo=weight{2};
% fclose(fid);
% 得到的是一个cell array,调用语句weight{1}、{2} ... , 以及weight{1}（1）。  
%
% TSinpo=interp1(XS,TS,Xinpo,'linear');     % interpolation with linear method
