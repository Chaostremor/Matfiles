% read data

% y = load('PGCopt_002_1.25-6.5Hz_2pass_40sps.input');
% ynorm = y/max(abs(y));
% sps = 40;
% x = linspace(0, floor((length(y)-1)/sps), length(y));
% figure
% plot(x, ynorm);

sps = 40;
y = opt(1:180*sps); 
ynorm = y/max(abs(y));
x = linspace(0, floor((length(y)-1)/sps), length(y));
figure (1)
plot(x, ynorm,'linewidth', 0.5, 'color', 'k');
xlabel('Time (s)', 'Fontsize', 12);
ylabel('Normailized amplitude', 'Fontsize', 12);
set(gca,'FontSize',10)
set(gca, 'xlim', [0, 180]);
set(gca, 'xtick', 0: 20: 180);
set(gca, 'ylim', [-1.1, 1.1]);
set(gca, 'ytick', -1: 0.5: 1.0);
scrsz=get(0,'ScreenSize');
set(gcf, 'unit', 'normalized', 'position', [0.05, 0.05, 0.3, 0.1], 'PaperPositionMode', 'auto');
% set(gcf, 'unit', 'centimeters', 'position', [5, 3, 20, 5], 'PaperPositionMode', 'auto');
print('-dpdf','-r600', [getenv('ALLAN'),'/waveexample.pdf']);
% print('-dpdf','-r600','-bestfit', [getenv('ALLAN'),'/waveexample.pdf']);