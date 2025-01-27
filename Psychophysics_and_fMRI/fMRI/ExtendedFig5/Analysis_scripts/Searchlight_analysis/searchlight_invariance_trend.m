%% A-P by plane analysis
% Many previous works show that there is increasing invariance from
% posterior to anterior regions of the HLVC (IT) cortex.

%To show this pattern, in-plane values are averaged and plotted from the
%anterior to posterior direction (Supplementary Fig. 5c)

clear
close all
clc

addpath('/isilon/LFMI/VMdrive/Software/r2020b/toolbox/images/iptformats/');

allsubs = {'001','004','005','006','007','008','009','010','012','013'};

filedir = '~/All_searchlight/LOandFG_searchlight_allsubs_3per_2mm.nii.gz';

subj_output = niftiread(filedir);
subj_output2 = double(subj_output);
maxval = max(max(max(max(subj_output2))));
for s = 1:10
    singlesub{s} = subj_output2(:,~:,:,s);
end

clearvars smoothZ

figure;
hold on
for s = 1:10
%     subplot(2,4,s)
    single_SLmap = singlesub{s};
    single_SLmap(single_SLmap==0)=NaN;
    allsub_SLmap(:,:,:,s) = single_SLmap;
    % Remove the first (X) and third (Z) dimensions by averaging:
    averageX = nanmean(single_SLmap,1);
    averageZ(s,:) = nanmean(averageX,3);
    for plane = 1:size(single_SLmap,2)
        maxdist(plane) = max(max(single_SLmap(:,plane,:)));
        meandist(plane) = nanmean(nanmean(single_SLmap(:,plane,:)));
    end
    meandist2 = meandist(~isnan(meandist));
    if length(meandist2) == 41
        meandist2 = [meandist2 meandist2(end)]; %add one value
    elseif length(meandist2) == 43
        meandist2 = meandist2(1:end-1); %remove last value if it's too long
    end
    for perm = 1:50
        tempperm(perm,:) = meandist2(randperm(length(meandist2)));
    end

    shuffledZ(s,:) = nanmean(tempperm);
    smoothZ(s,:) = movmean(meandist,6);
    smoothShuffledZ(s,:) = movmean(shuffledZ(s,:),6);
    maxZ(s,:) = maxdist;
    averageZ2(s,:) = meandist;

    plot(averageZ(s,:));
end

numSubjects = 10;

for s = 1:length(allsubs)
    sub = allsubs{s};
    for p = 1:1000 %500 # of perms
        clearvars perm_data
        perm_loc = fullfile(['~/All_searchlight/LOandFG/perm/',sub,'/',num2str(p-1),'_LOandFG_3per_2mm.nii.gz']);
        permdata = niftiread(perm_loc);
        permdata(permdata == 0) = NaN;
        averageX = nanmean(permdata,1);
        averageZ = nanmean(averageX,3);
        permdist = averageZ(~isnan(averageZ));
        permVals(s,p,:) = movmean(permdist,6);
    end
end

for r = 1:size(permdist,2) % for each plane
    for s = 1:numSubjects
        for p = 1:10000
            permvals_allsub = permVals(:,:,r);
            if ~isnan(nanmean(nanmean(permvals_allsub)))
                temp = randsample(permvals_allsub(:),numSubjects);
                perm_avgsubj(r,p,:) = mean(temp);
            else
                perm_avgsubj(r,p,:) = nan;
            end
        end
    end
end
clearvars permVals
for r = 1:size(permdist,2) % for each plane
    histperm = perm_avgsubj(r,:,:);
    sorted_histperm = sort(histperm);
    upperceiling_upper(r) = sorted_histperm(10000);
    upperceiling_med(r) = median(sorted_histperm);
    upperceiling_lower(r) = sorted_histperm(500);
end

% clearvars perm_avgsubj

%% t-test against 0
[h,pval_ttest] = ttest(smoothZ,0);
pval_ttest = pval_ttest(~isnan(pval_ttest));
corrected_pval_test = mafdr(pval_ttest,'BHFDR',true);

%% plotting
ap_vals = nanmean(smoothZ);
ap_vals_SEM = nanstd(smoothZ)/sqrt(10);
ap_vals_SHUFFLE = nanmean(smoothShuffledZ);
ap_vals_SEM_SHUFFLE = nanstd(smoothShuffledZ)/sqrt(10);

ap_vals = ap_vals(~isnan(ap_vals));
ap_vals_SEM = ap_vals_SEM(~isnan(ap_vals_SEM));

%ap_vals_SHUFFLE = ap_vals_SHUFFLE(~isnan(ap_vals_SHUFFLE));
%ap_vals_SEM_SHUFFLE = ap_vals_SEM_SHUFFLE(~isnan(ap_vals_SEM_SHUFFLE));

ap_vals = ap_vals(6:end-6);
ap_vals_SEM = ap_vals_SEM(6:end-6);
upperceiling_upper = upperceiling_upper(6:end-6);
upperceiling_med = upperceiling_med(6:end-6);
upperceiling_lower = upperceiling_lower(6:end-6);
% 
ap_vals_SHUFFLE = ap_vals_SHUFFLE(1:26);
ap_vals_SEM_SHUFFLE = ap_vals_SEM_SHUFFLE(1:26);
upperceiling_upper = upperceiling_upper(1:26);
upperceiling_med = upperceiling_med(1:26);
upperceiling_lower = upperceiling_lower(1:26);

perm_avgsubj2 = perm_avgsubj(6:end-6,:,:);
figure; set(gcf,'Color','w'); hold on
for r = 1:26
    subplot(4,7,r)
    histperm = perm_avgsubj2(r,:,:);
    sorted_histperm = sort(histperm);
    upperceiling_upper2 = sorted_histperm(10000);
    upperceiling_lower2 = sorted_histperm(500);

    pval = sum(histperm < ap_vals(:,r)) / numel(histperm);
    z_score(r) = (ap_vals(:, r) - mean(histperm)) / (std(histperm));
    histogram(histperm);
    xline(ap_vals(:,r),'r');
    xline(upperceiling_upper2);
    xline(upperceiling_lower2);
    if pval < 0.05
        title('invariant');
        invarplane(r) = 1;
        pvalall(r) = pval;
    else
        invarplane(r) = 0;
        pvalall(r) = pval;
    end
end

%% Cluster-based permutation test
thresh = 1.64;
[clust_size,start,stop,clust_labels] = find_clusters(abs(z_score) > thresh);

num_perm = 1000;
for p = 1:num_perm
    shuffled_binary = clust_labels(randperm(length(clust_labels)));
    [clust_size_shuffled,~,~,~] = find_clusters(shuffled_binary);
    permuted_max_cluster(p) = max(clust_size_shuffled);
end

sig_cluster = size(clust_size,1);
for c = 1:length(clust_size)
    clust_stats_pval(c) = sum(permuted_max_cluster > clust_size(c)) / num_perm;
    if clust_stats_pval<0.05
        sig_cluster(c) = 1; 
        sig_clust_start{c} = start(c);
        sig_clust_end{c} = stop(c);
    end
end

%% Plot parameters
greencolor = [0.4660 0.6740 0.1880];
shadeColor = [0.9290 0.6940 0.1250];  % RGB values for gray
alphaValue = 0.25;  % Transparency value (0 for fully transparent, 1 for fully opaque)

%% Plot with signifiance bar only:
figure; set(gcf,'Color','w');hold on
x = 1:length(ap_vals);
curve1 = ap_vals + ap_vals_SEM;
curve2 = ap_vals - ap_vals_SEM;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, greencolor);
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.2);
plot(ap_vals,'Color',greencolor,'LineWidth',1);
bold_indices = find(invarplane == 1);
% plot(bold_indices, ap_vals(bold_indices), 'Color', greencolor, 'LineWidth', 4);
set(gca,'XTick',[])
y_significance = 0.0001;  

for c = length(clust_size)
    start_idx =  sig_clust_start{c};
    end_idx =  sig_clust_end{c};
    plot([start_idx, end_idx], [y_significance, y_significance], 'Color',greencolor, 'LineWidth', 5);
end

%% Trend Analysis stats 1: Mann-Kendall
[H,p_value]=Mann_Kendall(ap_vals,0.05);

[H_shuffle,p_value_shuffle]=Mann_Kendall(ap_vals_SHUFFLE,0.05);

%% Trend Analysis stats 2: Fit a linear regression model
model = fitlm(1:numel(ap_vals), ap_vals);
disp(model)

model_shuffle = fitlm(1:numel(ap_vals_SHUFFLE), ap_vals_SHUFFLE);
disp(model_shuffle)

%% Plot with yellow ribbon:
%figure; set(gcf,'Color','w');
%hold on
%x = 1:length(ap_vals);
%curve1 = ap_vals + ap_vals_SEM;
%curve2 = ap_vals - ap_vals_SEM;
%x2 = [x, fliplr(x)];
%inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween, greencolor);
%set(h, 'EdgeColor', 'none'); % Remove black border
%alpha(0.3);
%plot(ap_vals,'Color',greencolor,'LineWidth',2);
%set(gca,'XTick',[])
%curve1 = ap_vals_SHUFFLE + ap_vals_SEM_SHUFFLE;
%curve2 = ap_vals_SHUFFLE - ap_vals_SEM_SHUFFLE;
%x2 = [x, fliplr(x)];
%inBetween = [curve1, fliplr(curve2)];
%h = fill(x2, inBetween,[0.7 0.7 0.7]);
%set(h, 'EdgeColor', 'none'); % Remove black border
%alpha(0.3);
%plot(ap_vals_SHUFFLE,'Color',[0.7 0.7 0.7],'LineWidth',2);

% figure; hold on
%roivals = 1:length(upperceiling_upper);
%X=[roivals,fliplr(roivals)]; %create continuous x value array for plotting
%Y=[upperceiling_lower,fliplr(upperceiling_upper)]; %create y values for out and then back
%fill(X,Y,shadeColor, 'FaceAlpha', alphaValue,'EdgeColor', 'none');   
%plot(roivals,upperceiling_med,'Color',shadeColor,'LineWidth',1.5);

