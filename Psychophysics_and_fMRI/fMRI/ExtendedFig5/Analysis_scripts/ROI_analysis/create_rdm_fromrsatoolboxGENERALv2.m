%% Across subject analysis
% Run rsa_plots_crossnobis_resamplerun.py first.
% This script will run permutations & plot the final figure (Supplementary Fig. 5)

clear all
close all
manipType = [{'Original'},{'CW'},{'CCW'},{'Rshift'},{'Lshift'},{'LRinvert'},{'Size'}];
distanceType = 3;
numManip = 7;

subjectList = [{'001'}, {'004'}, {'005'}, {'006'}, {'007'}, {'008'}, {'009'}, {'010'},{'012'},{'013'}];
numSubjects = length(subjectList);
roinames = {'new_v1_2highres' 'new_v2_2highres' 'new_v3_2highres' 'new_hv4_2highres'...
'new_LO1_2highres' 'new_LO2_2highres' 'FG_3percent_2highres'...
'new_IPS0_2highres' 'new_IPS1_2highres'...
'new_IPS2_2highres' 'new_IPS3_2highres' 'new_IPS4_2highres' 'new_IPS5_2highres' 'new_SPL1_2highres'...
'FPN_F_2highres' 'FPN_P_2highres' 'DMN_P_2highres' 'DMN_PCC_2highres' 'DMN_MPFC_2highres'};
roinames2 = {'V1' 'V2' 'V3' 'V4'...
'LO1' 'LO2' 'FC'...
'IPS0' 'IPS1'...
'IPS2' 'IPS3' 'IPS4' 'IPS5' 'SPL'...
'FPN-F' 'FPN-P' 'DMN-P' 'DMN-PCC' 'DMN-mPFC'};


grproi=input('Plotting grouped ROIs? [y/n]: ','s');

if grproi == 'y'
    subjectList = [{'001'},{'004'}, {'005'}, {'006'}, {'007'}, {'009'}, {'010'},{'012'},{'013'}];
    numSubjects = length(subjectList);
    roinames = {'EVC_2highres','OBJ_2highres','DORSAL_2highres','FPN_2highres','DMN_2highres'};
    %for figure labeling:
    roinames2 = {'EVC' 'OBJ' 'DORSAL' 'FPN' 'DMN'};
end

numROIs = length(roinames);


for s = 1:numSubjects
m = num2str(subjectList{s});
    for r = 1:length(roinames)
        roi = roinames{r};
        rdmALL{s,r} = load(['/Analysis/fMRI_analysis/Fig3_neural_dist_analysis/ROI_analysis_scripts/data/' m '/rdms/' m '_' roi '_crossnobisBYIMG_pseudo_rdm.mat'],'rdm');
        rdmbyROI {s,r} = reshape(rdmALL{s,r}.rdm,[numManip*10,numManip*10]); %rearrange by ROI, not by subject
    end
end

%Permutation testing per subject
ind = numManip;
imgbyimgMAT = cell(numROIs, 1);
for s = 1:numSubjects
    for r = 1:numROIs
        counter = 0;
        for i = 1:10
            for j = 1:10
                submatrix = rdmbyROI{s,r}((i-1)*ind+1:i*ind, (j-1)*ind+1:j*ind);
                if i == j
                    byimgALL_matrix{s,r}(:,:,i) = submatrix;
                elseif i > j
                    counter = counter+1;
                    submatrix = submatrix - diag(diag(submatrix)); %make diagonals zero to only take the lower triangle
                    %because naturally, submatrix-off-digonal matrices have non-zero
                    %diagonals.
                    temp = tril(submatrix);
                    permutedMatrix{s,r}(:,counter) = temp(temp~=0);
                end
            end
        end
    end
end

%average across all subjects, for every roi

for s = 1:numSubjects
    for r = 1:numROIs
        temp = byimgALL_matrix{s,r};
        averagedImg = mean(temp,3);
        temp2 = tril(averagedImg);
        temp3 = temp2(temp2~=0);
        subtrianglepersubj(s,r)=mean(temp3);
    end
end

for r = 1:numROIs
    subtriangle(r) = mean(subtrianglepersubj(:,r));
    [~,p(r)] = ttest(subtrianglepersubj(:,r),0);
    subtriangle_sem(r) = std(subtrianglepersubj(:,r))/sqrt(numSubjects);
end

correctedpval = mafdr(p,'BHFDR',true);
display(correctedpval);

figure;
hold on
bar(subtriangle);
plot(subtrianglepersubj','linestyle','none','marker','o','MarkerFaceColor', 'k');

%%

clearvars permDist perm2

for r = 1:numROIs
    temp = rdmbyROI(:,r);
    temp2 = cat(3,temp{:});
    rdmbyROI_avg{r} = mean(temp2,3);
end

for r = 1:numROIs
    for p = 1:10000
        for s = 1:numSubjects
            perm2 = cell2mat(permutedMatrix(s,r));
            for d = 1:(numManip*(numManip-1)/2)
                RANDimgselection = randsample(45,10);
                pairing(d,:) = perm2(d,RANDimgselection);
            end
            perm_all(p,s) = mean(pairing(:)); %10 images x 21 lower triangle vals
        end
    end
    
    for p = 1:10000
         ind = randperm(size(perm_all,1));
         ind = ind(1:numSubjects);
         temp = zeros(numSubjects,1);
         for s = 1:numSubjects
            temp(s) = perm_all(ind(s),s);
         end
         perm_avgsubj{r}(p,:) = mean(temp);
    end
end

figure;set(gcf,'color','w');
for r = 1:numROIs
    histperm = perm_avgsubj{r}(:);
    sorted_histperm = sort(histperm);
    
    upperceiling_upper(r) = sorted_histperm(10000);
    upperceiling_med(r) = median(sorted_histperm);
    upperceiling_lower(r) = sorted_histperm(500);
    
    subplot(4,7,r)
    hold on
    histogram(histperm);
    x_vertical = subtriangle(r);
    ylim = get(gca,'ylim');  % Get y-axis limits to determine line length
    line([x_vertical, x_vertical], ylim, 'Color', 'red');  % Draw vertical line
    
    line([upperceiling_upper(r), upperceiling_upper(r)], ylim, 'Color', 'black');  % Draw vertical line
    line([upperceiling_lower(r), upperceiling_lower(r)], ylim, 'Color', 'black');  % Draw vertical line
    
    pval = sum(histperm < x_vertical) / numel(histperm);
    p_value(r) = pval;
    title(roinames{r});
    xlabel(pval);
end

correctedpvalue = mafdr(p_value,'BHFDR',true);
display(correctedpvalue);

%% final fig plot
shadeColor = [0.9290 0.6940 0.1250];  % RGB values for gray
alphaValue = 0.25;  % Transparency value (0 for fully transparent, 1 for fully opaque)

clearvars ylim
figure;
set(gcf,'color','w');
hold on
p1 = patch([0 4.5 4.5 0], [max(upperceiling_upper + 0.001) max(upperceiling_upper+ 0.001) -0.002 -0.002], 'g','FaceAlpha',0.1,'EdgeColor', 'none');
p2 = patch([4.5 7.5 7.5 4.5], [max(upperceiling_upper + 0.001) max(upperceiling_upper+ 0.001) -0.002 -0.002], 'r','FaceAlpha',0.1,'EdgeColor', 'none');
p3 = patch([7.5 14.5 14.5 7.5], [max(upperceiling_upper+ 0.001) max(upperceiling_upper+ 0.001) -0.002 -0.002], 'b', 'FaceAlpha',0.1,'EdgeColor', 'none');
p4 = patch([14.5 16.5 16.5 14.5], [max(upperceiling_upper+ 0.001) max(upperceiling_upper+ 0.001) -0.002 -0.002], 'm', 'FaceAlpha',0.1,'EdgeColor', 'none');
p5 = patch([16.5 24.5 24.5 16.5], [max(upperceiling_upper+ 0.001) max(upperceiling_upper+ 0.001) -0.002 -0.002], 'c', 'FaceAlpha',0.1,'EdgeColor', 'none');

b = bar(subtriangle);
b.FaceColor = [0.4660 0.6740 0.1880];
b.EdgeColor = [0.4660 0.6740 0.1880];
h = errorbar(subtriangle,subtriangle_sem,'.');
h.Color = 'k'; h.LineWidth = 2;

roivals=1:numROIs;

X=[roivals,fliplr(roivals)]; %create continuous x value array for plotting
Y=[upperceiling_lower,fliplr(upperceiling_upper)]; %create y values for out and then back
fill(X,Y,shadeColor, 'FaceAlpha', alphaValue,'EdgeColor', 'none');   

plot(roivals,upperceiling_upper,'Color',shadeColor);
plot(roivals,upperceiling_med,'Color',shadeColor,'LineWidth',1.5);
plot(roivals,upperceiling_lower,'Color',shadeColor);

maxylimval=max(upperceiling_upper)+0.001;
xlim([0 20]);
ylim([-0.002 maxylimval]);
ylabel('Averaged lower triangle distances');

set(gca,'XTickLabel',roinames2, 'xtick',1:numel(roinames));
xtickangle(45);
ax = gca;
ax.FontSize = 15;
ylabel('Averaged Lower Triangle distances');
set(gca,'fontname','arial')  % Set it to times
hold off

pval = [0.244,0.1585,0.0583,0.0526,0.0103,0.0005,0.0097,0.0602,0.2272,0.251,0.3004,0.393,0.2962,0.5335,0.0122,0.2835,0.1922,0.8802,0.4417];
correctedpval = mafdr(pval,'BHFDR',true);
display(correctedpval);
