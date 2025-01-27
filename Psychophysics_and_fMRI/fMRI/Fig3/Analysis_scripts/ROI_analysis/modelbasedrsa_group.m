%% Plot group analysis for model-based RSA (Fig. 3):

clear all
close all
clc

cd 'Mooney_PriorInvariance/InvarianceMap_fmri/Data/model_compare_results/Spearman_H2d/'

grproi=input('Plotting grouped ROIs? [y/n]: ','s');
subjectList = [{'001'},{'004'}, {'005'}, {'006'}, {'007'}, {'008'},{'009'}, {'010'},{'012'},{'013'}];
if grproi == 'y'
    numSubjects = length(subjectList);
    roinames = {'EVC_2highres','OBJ_2highres','DORSAL_2highres','FPN_2highres','DMN_2highres'};
    %for figure labeling:
    roinames2 = {'EVC' 'OBJ' 'DORSAL' 'FPN' 'DMN'};
else
roinames = {'new_v1_2highres' 'new_v2_2highres' 'new_v3_2highres' 'new_hv4_2highres'...
        'new_LO1_2highres' 'new_LO2_2highres' 'FG_3percent_2highres'...
        'new_IPS0_2highres' 'new_IPS1_2highres'...
        'new_IPS2_2highres' 'new_IPS3_2highres' 'new_IPS4_2highres' 'new_IPS5_2highres' 'new_SPL1_2highres'...
        'FPN_F_2highres' 'FPN_P_2highres' 'DMN_P_2highres' 'DMN_PCC_2highres' 'DMN_MPFC_2highres'};
roinames2 = {'V1' 'V2' 'V3' 'V4'...
    'LO1' 'LO2' 'FG'...
    'IPS0' 'IPS1'...
    'IPS2' 'IPS3' 'IPS4' 'IPS5' 'SPL'...
    'FPN-F' 'FPN-P' 'DMN-P' 'DMN-PCC' 'DMN-mPFC'};
end

compare_results = cell(1, length(roinames));
for r = 1:length(roinames)
    data = load([roinames{r}, '_modelcomparison.mat']);
    compare_results{r} = data;
    mean_vals(r) = mean(data.evaluations);
    var_vals(r) = std(data.evaluations)/sqrt(10);
    lowerceiling(r) = data.noise_ceiling(1);
    upperceiling(r) = data.noise_ceiling(2);
end
%%
shadeColor = [0.9290 0.6940 0.1250];  % RGB values for gray
alphaValue = 0.25;  % Transparency value (0 for fully transparent, 1 for fully opaque)

figure;
set(gcf,'color','w');
hold on

p1 = patch([0 4.5 4.5 0], [max(upperceiling + 0.001) max(upperceiling+ 0.001) 0.0 0], 'g','FaceAlpha',0.1,'EdgeColor', 'none');
p2 = patch([4.5 7.5 7.5 4.5], [max(upperceiling + 0.001) max(upperceiling+ 0.001) 0 0], 'r','FaceAlpha',0.1,'EdgeColor', 'none');
p3 = patch([7.5 14.5 14.5 7.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) 0 0], 'b', 'FaceAlpha',0.1,'EdgeColor', 'none');
p4 = patch([14.5 16.5 16.5 14.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) 0 0], 'm', 'FaceAlpha',0.1,'EdgeColor', 'none');
p5 = patch([16.5 24.5 24.5 16.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) 0 0], 'c', 'FaceAlpha',0.1,'EdgeColor', 'none');

b = bar(mean_vals);
b.FaceColor = [0.4660 0.6740 0.1880];
b.EdgeColor = [0.4660 0.6740 0.1880];
h = errorbar(mean_vals,var_vals,'.');
h.Color = 'k'; h.LineWidth = 2;

roivals=1:length(roinames);

X=[roivals,fliplr(roivals)]; %create continuous x value array for plotting
Y=[lowerceiling,fliplr(upperceiling)]; %create y values for out and then back
fill(X,Y,shadeColor, 'FaceAlpha', alphaValue,'EdgeColor', 'none');

plot(roivals,upperceiling,'Color',shadeColor);
% plot(roivals,upperceiling_med,'Color',shadeColor,'LineWidth',1.5);
plot(roivals,lowerceiling,'Color',shadeColor);

maxylimval=max(upperceiling)+0.001;
xlim([0 25]);
% ylim([-0.03 0.01]);
ylabel('Spearmans rho');

set(gca,'XTickLabel',roinames2, 'xtick',1:numel(roinames));
xtickangle(45);
ax = gca;
ax.FontSize = 15;

hold off

x = 1:length(roinames); % X values for the bars
y = mean_vals; % Height of each bar

%% Create the bar plot
figure; set(gcf,'color','w');
hold on
p1 = patch([0 4.5 4.5 0], [max(upperceiling + 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'g','FaceAlpha',0.1,'EdgeColor', 'none');
p2 = patch([4.5 7.5 7.5 4.5], [max(upperceiling + 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'r','FaceAlpha',0.1,'EdgeColor', 'none');
p3 = patch([7.5 14.5 14.5 7.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'b', 'FaceAlpha',0.1,'EdgeColor', 'none');
p4 = patch([14.5 16.5 16.5 14.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'm', 'FaceAlpha',0.1,'EdgeColor', 'none');
p5 = patch([16.5 24.5 24.5 16.5], [max(upperceiling+ 0.001) max(upperceiling+ 0.001) -0.01 -0.01], 'c', 'FaceAlpha',0.1,'EdgeColor', 'none');

b = bar(x, y);
barWidth = b.BarWidth;

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

set(gca,'XTickLabel',roinames2, 'xtick',1:numel(roinames));
xtickangle(45);
ax = gca;
ax.FontSize = 15;

hold off; 

%% Evaluate perms
for r = 1:length(roinames)
    permdir = fullfile(['Mooney_PriorInvariance/InvarianceMap_fmri/Data/model_compare_results/Spearman_H2d/',roinames{r},'/perms']);
    for p = 1:100
        data = load(fullfile(permdir, [roinames{r}, '_modelcomparison_perm', num2str(p), '.mat']));
        perm_matrix{r}(p,:) = squeeze(data.evaluations);
    end
end

niterations = 50000;
nsubjects = length(subjectList);
figure;
for r = 1:length(roinames)
    matrix = perm_matrix{r};
    for pp = 1:niterations
        perm = randsample(100,10);
        for s = 1:nsubjects
            nullvals(s) = matrix(perm(s),s);
        end
        null(pp) = mean(nullvals);
    end
    subplot(4,7,r);
    histogram(null);
    ylim = get(gca,'ylim');  % Get y-axis limits to determine line length
    line([mean_vals(r), mean_vals(r)], ylim, 'Color', 'red');  % Draw vertical line
    title(roinames2(r));
    pvalues(r) = sum(mean_vals(r) <= null) / niterations;
    xlabel(pvalues(r));
end

correctedpvalue = mafdr(pvalues,'BHFDR',true);
display(correctedpvalue);
