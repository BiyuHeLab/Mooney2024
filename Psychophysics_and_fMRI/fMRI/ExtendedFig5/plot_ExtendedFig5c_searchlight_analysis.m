%% A-P by plane analysis (Supplementary Fig. 5c)
% Many previous works show that there is increasing invariance from
% posterior to anterior regions of the IT cortex.

%To show this pattern, slice the brain from the posterior to anterior axis, 
%take within-plane averages and plot the trend.

clear
close all
clc

addpath(genpath('~/.matlab/R2020b/vistasoft-master/'));
addpath('~/Software/r2020b/toolbox/images/iptformats/');

allsubs = {'001','004','005','006','007','008','009','010','012','013'};

filedir = '~/ExtendedFig5/searchlight_results_data/LOandFG_searchlight_allsubs_3per_2mm.nii.gz';

subj_output = niftiread(filedir);
subj_output2 = double(subj_output);
maxval = max(max(max(max(subj_output2))));
for s = 1:10
    singlesub{s} = subj_output2(:,:,:,s);
end

figure;
hold on
for s = 1:10
%     subplot(2,4,s)
    single_SLmap = singlesub{s};
    single_SLmap(single_SLmap==0)=NaN;
    allsub_SLmap(:,:,:,s) = single_SLmap;
    % Average across the first (X) and third (Z) dimensions:
    averageX = nanmean(single_SLmap,1);
    averageZ(s,:) = nanmean(averageX,3);
    for plane = 1:size(single_SLmap,2)
        maxdist(plane) = max(max(single_SLmap(:,plane,:)));
        meandist(plane) = nanmean(nanmean(single_SLmap(:,plane,:)));
    end

    meandist2 = meandist(~isnan(meandist));
    if length(meandist2) == 41
        %pad with 0 bc of slight voxel mismatch after transformation
        meandist2 = [meandist2 0 0]; 
    elseif length(meandist2) == 42
        meandist2 = [meandist2 0]; %pad with 0
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

%%% real data:
ap_vals = nanmean(smoothZ);
ap_vals_SEM = nanstd(smoothZ)/sqrt(10);
ap_vals_SHUFFLE = nanmean(smoothShuffledZ);
ap_vals_SEM_SHUFFLE = nanstd(smoothShuffledZ)/sqrt(10);

ap_vals = ap_vals(~isnan(ap_vals));
ap_vals_SEM = ap_vals_SEM(~isnan(ap_vals_SEM));

ap_vals_SHUFFLE = ap_vals_SHUFFLE(~isnan(ap_vals_SHUFFLE));
ap_vals_SEM_SHUFFLE = ap_vals_SEM_SHUFFLE(~isnan(ap_vals_SEM_SHUFFLE));

ap_vals = ap_vals(6:end-6);
ap_vals_SEM = ap_vals_SEM(6:end-6);

%%% spatially shuffle distance values to create shuffled data:
ap_vals_SHUFFLE = ap_vals_SHUFFLE(1:26);
ap_vals_SEM_SHUFFLE = ap_vals_SEM_SHUFFLE(1:26);

%%% plot:
greencolor = [0.4660 0.6740 0.1880];
figure; set(gcf,'Color','w');
hold on
x = 1:length(ap_vals);
curve1 = ap_vals + ap_vals_SEM;
curve2 = ap_vals - ap_vals_SEM;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween, greencolor);
set(h, 'EdgeColor', 'none'); 
alpha(0.3);
plot(ap_vals,'Color',greencolor,'LineWidth',2);
set(gca,'XTick',[])
curve1 = ap_vals_SHUFFLE + ap_vals_SEM_SHUFFLE;
curve2 = ap_vals_SHUFFLE - ap_vals_SEM_SHUFFLE;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
h = fill(x2, inBetween,[0.7 0.7 0.7]);
set(h, 'EdgeColor', 'none'); 
alpha(0.3);
plot(ap_vals_SHUFFLE,'Color',[0.7 0.7 0.7],'LineWidth',2);

for c = length(clust_size)
    start_idx =  16 %sig_clust_start{c}; %value from prev. cluster based perm test
    end_idx =  22 %sig_clust_end{c}; %value from prev. cluster based perm test
    plot([start_idx, end_idx], [0.0001, 0.0001], 'Color',greencolor, 'LineWidth', 5);
end

% Statistics using Mann-Kendall trend analysis:
[H,p_value]=Mann_Kendall(ap_vals,0.05);
disp(p_value);

[H,p_value_shuffle]=Mann_Kendall(ap_vals_SHUFFLE,0.05);
disp(p_value_shuffle);

% Fit a linear regression model (alternaive to Mann-Kendall, if you want to make parametric assumptions):
model = fitlm(1:numel(ap_vals), ap_vals);
disp(model)

model_shuffle = fitlm(1:numel(ap_vals_SHUFFLE), ap_vals_SHUFFLE);
disp(model_shuffle)
