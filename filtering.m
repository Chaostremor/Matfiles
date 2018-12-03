function [ dat ] = filtering (data, delta, lowf, highf)
%FILTERING Summary of this function goes here
%   Detailed explanation goes here

% this program is aimed to filter data using a 2nd order 2-pass butterworth
% filter. You can also use other filters by changing butter to others 

order = 4; % smaller order can handle much narrower freq band

fs = 1/delta;
nyq = fs/2;
[B,A] = butter(order,[lowf/nyq highf/nyq]);
dat = filtfilt(B,A,data);

end