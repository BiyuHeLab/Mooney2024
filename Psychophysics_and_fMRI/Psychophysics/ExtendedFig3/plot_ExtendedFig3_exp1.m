clear;
%cd to data dir
addpath('~\Violinplot-Matlab-master');
load('correct_exp1.mat'); %data

totalImages = 219; 
subjNum = 33;
totalManip = 10;
manipTypes = [{'Original (12 dva)'};{'Size-small (6dva)'};{'Size-big (24 dva)'};{'VF shift'};{'LR Inversion'};{'90 degree rot'};{'M-bias'};{'P-bias'};{'Masked GS'};{'Catch Images'}];
exp = 1;

%% Plotting grayscale recognition only:
% mean by subject
for gM = 1:10
    for subj = 1:size(correct.manip{1,1},2) %for each subject
        for trial_type = 1:3 %3 because preM/GS/postM
            meanbysubj{gM,trial_type}(:,subj) = nanmean(correct.manip{gM,trial_type}(:,subj));
        end
    end
end

% for bar plot
clearvars vioPts
individPts(:,1) = meanbysubj{1,2};
individPts(:,2) = meanbysubj{10,2};
individPts(:,3) = meanbysubj{2,2};
individPts(:,4) = meanbysubj{3,2};
individPts(:,5) = meanbysubj{5,2};
individPts(:,6) = meanbysubj{6,2};
individPts(:,7) = meanbysubj{4,2};
individPts([30,10,26],:) = [];
gsmean = nanmean(individPts);
gssem = nanstd(individPts)/sqrt(length(individPts));

manipTypes = [{'Original'};{'Catch Images'};{'Size Reduction'};{'Size Enlargement'};{'LR Inversion'};{'90 degree rot'};{'VF shift'}];
manipTypes2 = [{''}; manipTypes];
figure;
hold on
bar(gsmean,'w');
xticklabels(manipTypes2);
ylabel('% correct');

e = errorbar(1:7,gsmean,gssem,'.');
e.CapSize = 0;
e.LineWidth = 1;
e.Color = [0 0 0];
subjNum=30;
x = [ones(1,subjNum) 2*ones(1,subjNum) 3*ones(1,subjNum) 4*ones(1,subjNum) 5*ones(1,subjNum)...
    6*ones(1,subjNum) 7*ones(1,subjNum)]; %x data needed to overlay the swarmchart, below
individDataPts = individPts(:)';
s = swarmchart(x,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', 'k'); %,'*','LineWidth',2,'Color','k');
s.XJitter = 'rand';
s.XJitterWidth = 0.35;
s.MarkerFaceAlpha = .18;
ax=gca;ax.LineWidth=1.2;ax.FontSize = 20; 
% title('Grayscale recognition');
yticklabels([0 50 100]);

%%
for gM = 1:totalManip
for s = 1:subjNum
    for imgIdx = 1:totalImages
        GyesCorrect.manip{gM,1} = nan(totalImages,subjNum);
        GyesCorrect.manip{gM,2} = nan(totalImages,subjNum);
        GyesCorrect.manip{gM,3} = nan(totalImages,subjNum);
        
        GnoCorrect.manip{gM,1} = nan(totalImages,subjNum);
        GnoCorrect.manip{gM,2} = nan(totalImages,subjNum);
        GnoCorrect.manip{gM,3} = nan(totalImages,subjNum);
    end
end
end

for gM = 1:totalManip
    for s = 1:subjNum
        for imgIdx = 1:totalImages
                
                if correct.manip{gM,2}(imgIdx,s) == 1 %recognized the grayscale 
                    GyesCorrect.manip{gM,1}(imgIdx,s) = correct.manip{gM,1}(imgIdx,s);
                    GyesCorrect.manip{gM,3}(imgIdx,s) = correct.manip{gM,3}(imgIdx,s);
                else
                    GyesCorrect.manip{gM,1}(imgIdx,s) = NaN;
                    GyesCorrect.manip{gM,3}(imgIdx,s) = NaN;
                end

                if correct.manip{gM,2}(imgIdx,s) == 0 %did NOT recognize the grayscale
                    GnoCorrect.manip{gM,1}(imgIdx,s) = correct.manip{gM,1}(imgIdx,s);
                    GnoCorrect.manip{gM,3}(imgIdx,s) = correct.manip{gM,3}(imgIdx,s);
                else
                    GnoCorrect.manip{gM,1}(imgIdx,s) = NaN;
                    GnoCorrect.manip{gM,3}(imgIdx,s) = NaN;
                end
                
        end
    end
end

% Calclualte mean and S.E.M for each grayscale manipulation group
%Change: not mean, but just collect the %correct as they are and plot all
%individually.

% First - collecting all responses without averaging anything
%PRE
for gM = 1:totalManip %for each manipulation
    for s = 1:subjNum
        x = GyesCorrect.manip{gM,1}(:,s);
        y = GnoCorrect.manip{gM,1}(:,s);
        Gyes_bySubj_preALL {gM,s} = x;
        Gno_bySubj_preALL {gM,s} = y;
    end
end

%POST
for gM = 1:totalManip %for each manipulation
    for s = 1:subjNum
        x = GyesCorrect.manip{gM,3}(:,s);
        y = GnoCorrect.manip{gM,3}(:,s);
        Gyes_bySubj_postALL {gM,s} = x;
        Gno_bySubj_postALL {gM,s} = y;
    end
end

% Second - now we average
%PRE
for gM = 1:totalManip %for each manipulation
    for s = 1:subjNum
        x = GyesCorrect.manip{gM,1}(:,s);
        y = GnoCorrect.manip{gM,1}(:,s);
        Gyes_bySubj_pre (gM,s) = nanmean(x);
        Gyes_bySubj_preLength(gM,s) = length(x(~isnan(x)));
        Gno_bySubj_pre (gM,s) = nanmean(y);
        Gno_bySubj_preLength(gM,s) = length(y(~isnan(y)));
    end
end

%POST
for gM = 1:totalManip %for each manipulation
    for s = 1:subjNum
        x = GyesCorrect.manip{gM,3}(:,s);
        y = GnoCorrect.manip{gM,3}(:,s);
        Gyes_bySubj_post (gM,s) = nanmean(x);
        Gyes_bySubj_postLength(gM,s) = length(x(~isnan(x)));
        Gno_bySubj_post (gM,s) = nanmean(y);
        Gno_bySubj_postLength(gM,s) = length(y(~isnan(y)));
    end
end

Gyes_meanPre = nanmean(Gyes_bySubj_pre,2);
Gyes_meanPost = nanmean(Gyes_bySubj_post,2);

Gno_meanPre = nanmean(Gno_bySubj_pre,2);
Gno_meanPost = nanmean(Gno_bySubj_post,2);

Gyes_semPre = nanstd(Gyes_bySubj_pre,1,2)/sqrt(size(Gyes_bySubj_pre,1));
Gyes_semPost = nanstd(Gyes_bySubj_post,1,2)/sqrt(size(Gyes_bySubj_post,1));

Gno_semPre = nanstd(Gno_bySubj_pre,1,2)/sqrt(size(Gno_bySubj_pre,1));
Gno_semPost = nanstd(Gno_bySubj_post,1,2)/sqrt(size(Gno_bySubj_post,1));

%% Figures/plotting
origColor = [99/255 169/255 95/255];
catchColor = [250/255, 125/255, 0];

figure;
mooneyLabel = [{'Pre'}; repmat({'Pre'}, subjNum-1, 1); {'Post'}; repmat({'Post'}, subjNum-1, 1)];
maniplist=[1 10 2 3 5 6 4];
for g = 1:7
    gM = maniplist(g);
    clearvars vioPts
    vioPts(:,2) = Gyes_bySubj_pre(gM,:);
    vioPts(:,1) = Gyes_bySubj_post(gM,:);
    
    subplot(2,4,g);
    hold on;
    box off
    if gM == 1
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', origColor);
    elseif gM == 2 || gM == 3 || gM == 4 || gM == 5 || gM == 6
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
    title(string(manipTypes(g)),'FontSize', 20);
    ax = gca;
    ax.TitleHorizontalAlignment = 'center';
    box off
end


