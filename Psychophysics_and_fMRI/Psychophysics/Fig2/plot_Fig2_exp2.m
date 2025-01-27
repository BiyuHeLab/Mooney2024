clear
clc
close all

%cd to data dir
addpath('~\Violinplot-Matlab-master');
load('correct_exp2.mat');


%% Plotting - verbal response
manipTypes = [{'Original (12 DVA)'};{'Mbias'};{'Pbias'};{'Catch'};{'Line drawing'};{'Category'}]; %for sake of titles
totalManip = 6;
subjNum = 12;

%mean by subject
for gM = 1:totalManip
    for subj = 1:size(correct.manip{1,1},2) %for each subject
        for trial_type = 1:3 %3 because preM/GS/postM
            meanbysubj{gM,trial_type}(:,subj) = nanmean(correct.manip{gM,trial_type}(:,subj));
        end
    end
end

%all subject mean
for gM = 1:totalManip
    for trial_type = 1:3 %3 because preM/GS/postM
        subjAllmean{gM,trial_type} = nanmean(meanbysubj{gM,trial_type});
    end
end

%% Plots for presentation/paper:
clearvars allmean2 allstd2 individDataPts meanPre meanPost semPre semPost
linedrawingcolor = [126/255 176/255 222/255];
origColor = [99/255 169/255 95/255];
catchColor = [250/255, 125/255, 0];

xbar = [1 2 4 5 7 8 10 11 13 14 16 17];
plotThis = [1 5 2 3 6 4];
manipTypes = [{'Original (12 DVA)'};{'Mbias'};{'Pbias'};{'Line drawing'};{'Category'};{'Catch'};];
x = [ones(1,subjNum) 2*ones(1,subjNum) 4*ones(1,subjNum) 5*ones(1,subjNum)...
    7*ones(1,subjNum) 8*ones(1,subjNum) 10*ones(1,subjNum) 11*ones(1,subjNum)...
    13*ones(1,subjNum) 14*ones(1,subjNum) 15*ones(1,subjNum) 16*ones(1,subjNum)]; %x data needed to overlay the swarmchart, below

% violin plot
mooneyLabel = [{'Pre'}; repmat({'Pre'}, 11, 1); {'Post'}; repmat({'Post'}, 11, 1)];
for g = 1:length(plotThis)
    clearvars vioPts
    gM = plotThis(g);
    vioPts(:,2) = meanbysubj{gM,1};
    vioPts(:,1) = meanbysubj{gM,3};

    % Jitter duplicate values: add rand()*0.1 when needed, but don't go
% above max val. If above max val, subtract rand()*0.1
[~, ~, uniqueval] = unique(vioPts, 'rows', 'stable');  
counts = accumarray(uniqueval, 1); %find how many times the unique val happens
repeat_ind = find(counts(uniqueval) > 1); % find indices of all repeat rows

    subplot(2,4,g)
    hold on;
    if gM == 1
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', origColor);
    elseif gM == 4
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', [250/255, 125/255, 0;250/255, 125/255, 0]);
    elseif gM == 3
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', linedrawingcolor);
    else
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', [1, 0.5, 0.7961;1, 0.5, 0.7961]);
    end
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
    fname = ({'Pre','Post'});
    set(gca, 'XTick', 1:length(fname),'XTickLabel',fname,'FontSize',16);
    ylim([0 1]);
    xlim([0.5 2.5]);
    ylabel ('% correct');
    title(string(manipTypes(g)),'FontSize', 20);
    ax = gca;
    ax.TitleHorizontalAlignment = 'center';
    box off
end
