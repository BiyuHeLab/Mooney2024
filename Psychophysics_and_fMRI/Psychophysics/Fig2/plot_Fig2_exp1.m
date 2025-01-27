clear;
clc
close all

%cd to data dir
addpath('~\Violinplot-Matlab-master');
load('correct_exp1.mat'); %data

%% Plotting 
manipTypes = [{'Original'};{'Size Reduction'};{'Size Enlargement'};{'VF shift'};{'LR Inversion'};{'90 degree rot'};{'M-bias'};{'P-bias'};{'Masked GS'};{'Catch Images'}];

%mean by subject
for gM = 1:10
    for subj = 1:size(correct.manip{1,1},2) %for each subject
        for trial_type = 1:3 %3 because preM/GS/postM
            meanbysubj{gM,trial_type}(:,subj) = nanmean(correct.manip{gM,trial_type}(:,subj));
        end
    end
end

%% violin plot

origColor = [99/255 169/255 95/255];
catchColor = [250/255, 125/255, 0];

plotThis = [1 2 3 4 5 6 10];
mooneyLabel = [{'Pre'}; repmat({'Pre'}, 29, 1); {'Post'}; repmat({'Post'}, 29, 1)];
for g = 1:length(plotThis)
    clearvars vioPts
    gM = plotThis(g);
    vioPts(:,2) = meanbysubj{gM,1};
    vioPts(:,1) = meanbysubj{gM,3};
    vioPts(3,:) = NaN;
    vioPts(10,:) = NaN;
    vioPts(26,:) = NaN;
    %remove nan rows:
    vioPts([30,10,26],:) = [];
    subplot(2,4,g);
    hold on;
    box off
    if g == 1
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', origColor);
    elseif g == 2 || g == 3 || g == 4 || g == 5 || g ==6
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', [1, 0.5, 0.7961;1, 0.5, 0.7961]);
    else
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', catchColor);
    end
    fname = ({'Pre','Post'});
    set(gca, 'XTick', 1:length(fname),'XTickLabel',fname,'FontSize',16);

    % Jitter duplicate values: add rand()*0.1 when needed, but don't go
    % above max val. If above max val, subtract rand()*0.1
    [~, ~, uniqueval] = unique(vioPts, 'rows', 'stable');  
    counts = accumarray(uniqueval, 1); %find how many times the unique val happens
    repeat_ind = find(counts(uniqueval) > 1); % find indices of all repeat rows
    
    maxval1 = max(vioPts(:,1)); minval1 = min(vioPts(:,1));
    maxval2 = max(vioPts(:,2)); minval2 = min(vioPts(:,2));
    maxval = min(maxval1, maxval2);
    minval = min(minval1, minval2);
    for i = 1:length(vioPts)
        if ismember(i, repeat_ind)
            temp = vioPts(i,:);
            epsilon = rand()*0.05;
            if temp(1) < maxval && temp(2) < maxval && temp(1) > minval && temp(2) > minval || temp(1)==0 || temp(2)==0
                vioPts(i,1) = temp(1) + epsilon;
                vioPts(i,2) = temp(2) + epsilon;
            else
                vioPts(i,1) = temp(1) - epsilon;
                vioPts(i,2) = temp(2) - epsilon;
            end            
        end
    end
    % Plot the modified values here:
    for i = 1:length(vioPts)
        plot([1, 2], [vioPts(i,2), vioPts(i,1)], '-', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.01);
    end
    ylim([0 1]);
    xlim([0.5 2.5]);
    ylabel ('% correct');
    title(string(manipTypes(gM)),'FontSize', 20);
    ax = gca;
    ax.TitleHorizontalAlignment = 'center';
    box off
end
