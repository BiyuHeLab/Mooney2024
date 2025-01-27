%% Mooney project analysis

%File labeling convention: all files are 4 digits long.
%The first number indicates file type (e.g. 0 for Mooney, 1 for Grayscale)
%Last three numbers is the image number (1 through 213).
%For example: 0041.bmp means Mooney ("0"), 41st image;
%1042.bmp means Grayscale ("1"), 42nd image.

clear;
close all;

%% Step 1: Import csv file; extract image file type, reaction time, response correctness

%cd to datapath
addpath('~\Violinplot-Matlab-master');
load('correct_exp2.mat');

%%%%%%%%%%%%%%%%%%%% Excel files from PsychoPy:
subjNum = 12; 

s = subjNum;
totalImages = 219;

% T = readtable('1c_mooney_control_August5_2022_2022_Aug_05_1543b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj1c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('2c_mooney_control_August5_2022_2022_Aug_08_1356b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj2c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('3c_mooney_control_August5_2022_2022_Aug_11_1210b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj3c_Resp.xlsx','PreserveVariableNames',true);
% % 
% T = readtable('4c_mooney_control_August5_2022_2022_Aug_12_1358b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj4c_Resp.xlsx','PreserveVariableNames',true);
% % 
% T = readtable('5c_mooney_control_August5_2022_2022_Aug_16_0958b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj5c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('6c_mooney_control_August5_2022_2022_Aug_19_1539b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj6c_Resp.xlsx','PreserveVariableNames',true);
% 
% T = readtable('7c_mooney_control_August5_2022_2022_Aug_24_1820b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj7c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('8c_mooney_control_August5_2022_2022_Aug_25_1341b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj8c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('9c_mooney_control_August5_2022_2022_Sep_06_1954b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj9c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('10c_mooney_control_August5_2022_2022_Oct_27_1101b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj10c_Resp.xlsx','PreserveVariableNames',true);

% T = readtable('11c_mooney_control_August5_2022_2022_Nov_03_0945b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
% T2 = readtable('Subj11c_Resp.xlsx','PreserveVariableNames',true);

T = readtable('12c_mooney_control_August5_2022_2022_Nov_11_1309b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
T2 = readtable('Subj12c_Resp.xlsx','PreserveVariableNames',true);

R3 = table2array(T2(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Images list (R1)

R1 = table2array(T(:,6)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make row number uniform, depending on valid blocks.
%300 per valid block.
matrixHeight = 900; %This is how many rows you should have, before you take away em

R1 = R1(1:matrixHeight);
T = T(1:matrixHeight,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Images ORDER list (R2)

%Issue -- sometimes, order is in 56th or 59 col, depending on trial
%Solution : Find 'order', or the last non-empty element, and take 6 rows below that.

R2 = cell(matrixHeight,1);
for ii = 1:matrixHeight
    if ismember(R1(ii,:),'img')
        lastCol = find(~cellfun(@isempty,table2array(T(ii,:))));
        R2(ii:ii+6,:) = table2array(T(ii:ii+6,lastCol(end)));
    end
end

R2 = str2double(R2(:))+1; %add 1's to order numbers so it's 1, 2, 3 not 0, 1, 2
R2 = R2(1:matrixHeight);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Rearrange images list (R1) based on order (R2)

%imgIdx is the CORRECT image list, based on presentation order in task.
imgIdx = zeros(length(R1),2);

for rr = 1:length(R1)
    %Step 1 - Extract image index number
    if contains(R1{rr},'img')
        imgIdx(rr,1) = 999; %999 will signify "img" or beginning of img chunk
    elseif contains(R1{rr},'bmp')
        str = R1{rr};
        fileNum = erase(str,'.bmp');
        check_img = fileNum(end-2:end); %Take the last 3 numbers, as image file number
        imgIdx(rr,1) = str2double(regexprep(check_img,'^0*','')); %take off leading zeros
        %Step 2 - Add on manipulation types next to image index column
            if length(fileNum) == 5
                imgIdx(rr,2) = 10;
            else
                imgIdx(rr,2) = str2double(fileNum(1));
            end
    else
        continue;
    end
end

%Sort those 2 columns based on "order", or R2

%Make NaNs in R2 zero
R2(isnan(R2))=0;

imgIdx2 = imgIdx(:,2); %just the manipulation types

sortedImgs = cell(length(imgIdx),2);

for ii = 1:length(imgIdx)-1
    if imgIdx(ii,1) == 999 %If the row is the start of an 'img' block
        tempImg = imgIdx(ii+1:ii+6)'; %collect img numbers here
        tempManip = imgIdx2(ii+1:ii+6)'; %collect manipulation types here
        tempOrder = R2(ii+1:ii+6); %collect img presentation order here
        
        %get rid of zeros
        tempImg = tempImg(tempImg~=0);
        tempOrder = tempOrder(tempOrder~=0);

        [~,OrderSort] = sort(tempOrder); %I think this was redundant, tempOrder alone should've been fine...
        sortedImgs{ii,1} = tempImg(OrderSort); %re-order img order based on presentation
        sortedImgs{ii,2} = tempManip(OrderSort); %re-order corresponding manip.
        
    end
end

%remove empty's
sortedImgs = sortedImgs(~cellfun(@isempty, sortedImgs(:,1)), :);
sortedImgsList = cell2mat(sortedImgs(:,1)); %turn images into list
sortedManipList = cell2mat(sortedImgs(:,2)');

sortedAll = [sortedImgsList sortedManipList'];

sortedAll(:,3)= R3; %This is correct/incorrect answers.

%Now sort by image number:
sortedAll2 = sortrows(sortedAll,1); %based on column 1 only

%Summary for "sortedAll":
%Col 1 - image index
%Col 2 - manipulation type
%Col 3 - correct/incorrect

%% Step 3: Analyze % correct for pre/GS/post, sorted by manipulation type

% Add correctness in column 2
totalImages = 219; %not 219 images, but image index labels range from 1-219
totalManip = 6;

mooney_yesno = repmat([0; 1; 0],matrixHeight/10,1); %0 if it's a mooney, 1 if it's a GS
mooney_prepost = repmat([1;0;0],matrixHeight/10,1); %1 if it's a pre-M, 0 if it's a post-M


for gM = 1:totalManip
    for rr = 1:totalImages
        correct.manip{gM,1}(:,subjNum) = nan(totalImages,1); %pre
        correct.manip{gM,2}(:,subjNum) = nan(totalImages,1); %grayscale
        correct.manip{gM,3}(:,subjNum) = nan(totalImages,1); %post
    end
end

for imgIdx = 1:length(sortedAll)
    if ~isnan(sortedAll2(imgIdx,1))
    if mooney_yesno(imgIdx) == 0 %if it's a Mooney image
        img = sortedAll2(imgIdx,1);
        if mooney_prepost(imgIdx) == 1 %1 means pre-Mooney
            gM = sortedAll2(imgIdx+1,2); %take the manipulation from the GS, 2nd col, one row down
            correct.manip{gM,1}(img,s) = sortedAll2(imgIdx,3); %3rd column is verbal response
        else %mooney_prepost is 0 means post-Mooney
            gM = sortedAll2(imgIdx-1,2); %take the manipulation from the GS, 2nd col, one row up
            correct.manip{gM,3}(img,s) = sortedAll2(imgIdx,3);
        end
    else
        gM = sortedAll2(imgIdx,2); %manipulation is the corresponding 2nd col
        correct.manip{gM,2}(img,s) = sortedAll2(imgIdx,3);
    end
    end
end

%% Important: save correct.mat

save('correct','correct');

%% Plotting - verbal response
manipTypes = [{'Original (12 DVA)'};{'Mbias'};{'Pbias'};{'Catch'};{'Line drawing'};{'Category'}]; %for sake of titles

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
catchColor = [127/255 0 255/255];

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


