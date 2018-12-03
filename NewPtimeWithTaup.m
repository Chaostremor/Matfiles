% test tauptime
clear;
dis = 30: 0.1: 97;
sdis = length(dis);
dep = 0: 1: 100;
sdep = length(dep);

%% way 1
% Ptt = load('C:\Users\Song Chao\Documents\MATLAB\myself\Ptime');
% Ptt = Ptt' ;
% save('C:\Users\Song Chao\Documents\MATLAB\myself\libBP\Ptimesdepthto100.mat', 'Ptt', 'dep', 'dis');


%% way 2
Ptt = zeros(sdep, sdis);

for idep = 1: sdep
    for idis = 1: sdis
        TT = tauptime('mod', 'prem', 'dep', dep(idep), 'ph', 'P', 'degrees', dis(idis));
        Ptt(idep, idis) = TT.time;
    end
end
save('C:\Users\Song Chao\Documents\MATLAB\myself\libBP\Ptimesdepthto100Matlab.mat', 'Ptt', 'dep', 'dis');

