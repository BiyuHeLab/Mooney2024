%% Plot Fig. 3c
clear; clc;

cd '~/Desktop/Ayaka/Mooney_PriorInvariance/InvarianceMap_fmri/Plot_sourcedata/Fig3/'

load('upperceiling.mat');
load('lowerceiling.mat');
load('mean_vals.mat');
load('var_vals.mat');
roinames = {'V1' 'V2' 'V3' 'V4'...
    'LO1' 'LO2' 'FG'...
    'IPS0' 'IPS1'...
    'IPS2' 'IPS3' 'IPS4' 'IPS5' 'SPL'...
    'FPN-F' 'FPN-P' 'DMN-P' 'DMN-PCC' 'DMN-mPFC'};

figure; set(gcf,'color','w');
hold on
p1 = patch([0 4.5 4.5 0], [max(upperceiling + 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'g','FaceAlpha',0.1,'EdgeColor', 'none');
p2 = patch([4.5 7.5 7.5 4.5], [max(upperceiling + 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'r','FaceAlpha',0.1,'EdgeColor', 'none');
p3 = patch([7.5 14.5 14.5 7.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'b', 'FaceAlpha',0.1,'EdgeColor', 'none');
p4 = patch([14.5 16.5 16.5 14.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'm', 'FaceAlpha',0.1,'EdgeColor', 'none');
p5 = patch([16.5 24.5 24.5 16.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'c', 'FaceAlpha',0.1,'EdgeColor', 'none');

x = 1:length(roinames); % X values for the bars
y = mean_vals; % Height of each bar

b = bar(x, y);
barWidth = b.BarWidth;

% Min and max y values for the rectangles
yMin = lowerceiling;
yMax = upperceiling;

hold on;

b = bar(x, y);
b.FaceColor = [0.4660 0.6740 0.1880];
b.EdgeColor = [0.4660 0.6740 0.1880];
h = errorbar(mean_vals,var_vals,'.');
h.Color = 'k'; h.LineWidth = 2;

roivals=1:length(roinames);
xlim([0 24]);
ylim([-0.01 0.017]);
ylabel('Spearmans rho');

set(gca,'XTickLabel',roinames, 'xtick',1:numel(roinames));
xtickangle(45);
ax = gca;
ax.FontSize = 15;
