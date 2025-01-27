%% Average & Plot:
load('meanbysubjVF.mat');

figure;
for gM = 1:3
    clearvars vioPts
    vioPts(:,2) = meanbysubjVF{gM,1};
    vioPts(:,1) = meanbysubjVF{gM,3};
    vioPts(3,:) = NaN;
    vioPts(10,:) = NaN;
    vioPts(26,:) = NaN;
    %remove nan rows:
    vioPts([30,10,26],:) = [];
    subplot(2,4,gM);
    hold on;
    box off
    if gM == 1
        vs = violinplot(vioPts(:), mooneyLabel,'ViolinColor', origColor);
    elseif gM == 2
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
    title(manipTypes(gM),'FontSize', 20);
    ax = gca;
    ax.TitleHorizontalAlignment = 'center';
    box off
end

%% statistics using generalized linear mixed model (GLMM) w/ logistic linking function

for m = 1:3
    for trialtype = [1 3]
        matrix = correctVFtrial(m,trialtype);
        for s = 1:size(matrix,2)
            for t = 1:size(matrix,1)
                temp = matrix{s,t};
                for r = 1:size(temp,2)
                    A = temp(:,r);
                    A = A(~isnan(A));
                    temp3{r} = A;
                end
            end
        end
        nonanmatrix{m,trialtype}=temp3;
    end
end

lmmpairs = {[1 2], [3 2]};
clc
for p = 1:length(lmmpairs)
    pairs = lmmpairs{p};
    lmmData = table();
    allrecograte = []; alltrialdata = []; allsubjdata = []; alltrialtype = [];
    allmaniptype = []; 
    for m = pairs
        % mcounter = mcounter + 1;
        for s = 1:length(nonanmatrix{1,1}) %for subject number 1 to 33
            for trialtype = [1 3]
                temp = nonanmatrix{m,trialtype}{s};
                allrecograte = [allrecograte; temp];
                trialdata = 1:length(temp);
                subjectdata = repmat(s,length(temp),1);
                alltrialdata = [alltrialdata; trialdata'];
                allsubjdata = [allsubjdata; subjectdata];
                if trialtype == 1
                    prepost = repmat({'pre'},length(temp),1);
                else
                    prepost = repmat({'post'},length(temp),1);
                end
                alltrialtype = [alltrialtype; prepost];
                manip = repmat(manipTypes(m),length(temp),1);
                allmaniptype = [allmaniptype; manip];
            end
        end
        lmmData = [lmmData;table(allsubjdata,alltrialdata,alltrialtype,allmaniptype,allrecograte,'VariableNames', {'subject', 'trialnum','trialtype', 'manip', 'recog_rate'})];
        lmm = fitglme(lmmData, 'recog_rate ~ manip * trialtype + (1|subject) + (1|trialnum)','Distribution','Binomial', 'Link','logit');
        if (length(lmm.Coefficients.pValue)) == 4
            pval_list2{p,1} = lmm.Coefficients.Name(end);
            pval_list2{p,2} =lmm.Coefficients.pValue(end);
        end
        disp(lmm);
    end
end


