clear;
close all;

%% Step 1: Import csv file; extract image file type, reaction time, response correctness

addpath('~\Violinplot-Matlab-master');
load('correct_3dvaVF.mat');

eyetrackcorrect = 'y';

%Add edf2mat path:
addpath('~\edf-converter-master\@Edf2Mat');

%%%%%%%%%%%%%%%%%%%% Excel files from PsychoPy:

totalImages = 219;    
    
%%everything else%%

for s = 1:33
    subjNum = s;
    if s == 10 || s == 26 || s == 3 %Exclude these participants
        for gM = 1:10
            for rr = 1:totalImages
                correct.manip{gM,1}(:,subjNum) = nan(totalImages,1); %pre
                correct.manip{gM,2}(:,subjNum) = nan(totalImages,1); %grayscale
                correct.manip{gM,3}(:,subjNum) = nan(totalImages,1); %post
            end
        end
    else
    if s == 1
        Ta = readtable('1_mooney_Main_March_23_2022_2022_Mar_28_1148b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1
        Tb = readtable('1b_mooney_Main_March_23_2022_2022_Mar_28_1221b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 1b
        T2 = readtable('Subj1_Resp.xlsx','PreserveVariableNames',true);
        edffile = '1b_001.edf';
    elseif s == 2
        T = readtable('2b_mooney_Main_March_23_2022_2022_Mar_29_1859b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 2
        T2 = readtable('Subj2_Resp.xlsx');
        edffile = '2b_001.edf';
    % elseif s == 3
    %     T = readtable('3_mooney_Main_March_23_2022_2022_Apr_01_1703b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 3
    %     T2 = readtable('Subj3_Resp.xlsx');
    % Exclude this participant who was clearly not paying attention, saying yes
    % to everything without verbal response.
    elseif s == 4
        T = readtable('4_mooney_Main_March_23_2022_2022_Apr_04_1909b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 4
        T2 = readtable('Subj4_Resp.xlsx');
        edffile = '4_001.edf';
    elseif s == 5
        T = readtable('5_mooney_Main_March_23_2022_2022_Apr_12_1418b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 5
        T2 = readtable('Subj5_Resp.xlsx');
        edffile = '5_002.edf';
    elseif s == 6
        T = readtable('6_mooney_Main_March_23_2022_2022_Apr_18_1931b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 6
        T2 = readtable('Subj6_Resp.xlsx');
        edffile = '6_001.edf';
    elseif s == 7
        T = readtable('7_mooney_Main_March_23_2022_2022_Apr_21_0858b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 7
        T2 = readtable('Subj7_Resp.xlsx');
        edffile = '7_001.edf';
    elseif s == 8
        T = readtable('8_mooney_Main_March_23_2022_2022_May_09_1256b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 8
        T2 = readtable('Subj8_Resp.xlsx');
        edffile = '8_001.edf';
    elseif s == 9
        T = readtable('9_mooney_Main_March_23_2022_2022_May_13_1024b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 9
        T2 = readtable('Subj9_Resp.xlsx');
        edffile = '9_002.edf';

    % %Subj 10 --> Bad conditions file, manipulations are cross matched (like
    % enlarged M-bias image :/) Unfortunatley, you cannot use this.
    elseif s == 11
        T = readtable('11_mooney_Main_March_23_2022_2022_May_16_1518b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 11
        T2 = readtable('Subj11_Resp.xlsx');
        edffile = '11_001.edf';
    elseif s == 12
        T = readtable('12_mooney_Main_March_23_2022_2022_May_18_1412b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 12
        T2 = readtable('Subj12_Resp.xlsx');
        edffile = '12_001.edf';
    elseif s == 13
        T = readtable('13_mooney_Main_March_23_2022_2022_May_31_1329b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 13
        T2 = readtable('Subj13_Resp.xlsx','PreserveVariableNames',true);
        edffile = '13_002.edf';
    elseif s == 14
        T = readtable('14_mooney_Main_March_23_2022_2022_Jun_06_1712b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 14
        T2 = readtable('Subj14_Resp.xlsx','PreserveVariableNames',true);
        edffile = '14_002.edf';
    elseif s == 15
        T = readtable('15_mooney_Main_March_23_2022_2022_Jun_07_1525b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 15
        T2 = readtable('Subj15_Resp.xlsx','PreserveVariableNames',true);
        edffile = '15_001.edf';
    elseif s == 16
        T = readtable('16_mooney_Main_March_23_2022_2022_Jun_08_1431b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 16
        T2 = readtable('Subj16_Resp.xlsx','PreserveVariableNames',true);
        edffile = '16_001.edf';
    elseif s == 17
        T = readtable('17_mooney_Main_March_23_2022_2022_Jun_08_1729b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 17
        T2 = readtable('Subj17_Resp.xlsx','PreserveVariableNames',true);
        edffile = '17_004.edf';
    elseif s == 18
        T = readtable('18_mooney_Main_March_23_2022_2022_Jun_09_1057b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 18
        T2 = readtable('Subj18_Resp.xlsx','PreserveVariableNames',true);
        edffile = '18_001.edf';
    elseif s == 19
        T = readtable('19_mooney_Main_March_23_2022_2022_Jun_15_1059b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 19
        T2 = readtable('Subj19_Resp.xlsx','PreserveVariableNames',true);
        edffile = '19_002.edf';
    elseif s == 20
        T = readtable('20_mooney_Main_March_23_2022_2022_Jun_20_1600b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 20
        T2 = readtable('Subj20_Resp.xlsx','PreserveVariableNames',true);
        edffile = '20_001.edf';
    elseif s == 21
        T = readtable('21_mooney_Main_March_23_2022_2022_Jun_23_1454b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 21
        T2 = readtable('Subj21_Resp.xlsx','PreserveVariableNames',true);
        edffile = '21_001.edf';
    elseif s == 22
        T = readtable('22_mooney_Main_March_23_2022_2022_Jun_27_1104b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 22
        % no eyetracking data
        T2 = readtable('Subj22_Resp.xlsx','PreserveVariableNames',true);
    elseif s == 23
        T =readtable('23_mooney_Main_March_23_2022_2022_Jun_30_1003b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 23
        T2 = readtable('Subj23_Resp.xlsx','PreserveVariableNames',true);
        edffile = '23_001.edf';
    elseif s == 24
        T = readtable('24_mooney_Main_March_23_2022_2022_Jul_01_1201b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 24
        T2 = readtable('Subj24_Resp.xlsx','PreserveVariableNames',true);
        edffile = '24_001.edf';
    elseif s == 25
        T = readtable('25_mooney_Main_March_23_2022_2022_Jul_06_1213b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 25
        T2 = readtable('Subj25_Resp.xlsx','PreserveVariableNames',true);
        edffile = '25_001.edf';
% Subj 26 - did not run task, because he said he was too tired.
    elseif s == 27
        T = readtable('27_mooney_Main_March_23_2022_2022_Jul_07_1002b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 27
        T2 = readtable('Subj27_Resp.xlsx','PreserveVariableNames',true);
        edffile = '27_001.edf';
    elseif s == 28
        T = readtable('28_mooney_Main_March_23_2022_2022_Jul_07_1249b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 28
        T2 = readtable('Subj28_Resp.xlsx','PreserveVariableNames',true);
        edffile = '28_001.edf';
    elseif s == 29
        T = readtable('29_mooney_Main_March_23_2022_2022_Jul_08_1122b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 29
        T2 = readtable('Subj29_Resp.xlsx','PreserveVariableNames',true);
        edffile = '29_001.edf';
    elseif s == 30
        T = readtable('30_mooney_Main_March_23_2022_2022_Jul_11_1640b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 30
        T2 = readtable('Subj30_Resp.xlsx','PreserveVariableNames',true);
        edffile = '30_001.edf';
    elseif s == 31
        T = readtable('31_mooney_Main_March_23_2022_2022_Jul_15_1349b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 30
        T2 = readtable('Subj31_Resp.xlsx','PreserveVariableNames',true);
        edffile = '31_001.edf';
    elseif s == 32
        T = readtable('32_mooney_Main_March_23_2022_2022_Jul_15_1538b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 32
        T2 = readtable('Subj32_Resp.xlsx','PreserveVariableNames',true);
        edffile = '32_001.edf';
    elseif s == 33
        T = readtable('33_mooney_Main_March_23_2022_2022_Jul_18_1000b1_loop.csv','ReadVariableNames',false,'PreserveVariableNames',true); %Subj 33
        T2 = readtable('Subj33_Resp.xlsx','PreserveVariableNames',true);
        edffile = '33_001.edf';
    end

%%%%%%%%%%%%%%%%%%%%% Unique modifications for certain subjects.

if s == 1
    %For subject 1, just combine T and Tb into one "T"
    Tb(1:height(Ta),:) = Ta;
    T = Tb;
end

%%%%%%%%%%%%%%%%%%%%% Subjective Responses (R3):

% T2 = readtable('Subj12_Resp.xlsx','PreserveVariableNames',true);
R3 = table2array(T2(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Images list (R1)

R1 = table2array(T(:,6)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make row number uniform, depending on valid blocks.
%300 per valid block.
matrixHeight = 900; %This is how many rows you should have, before you take away em
if s == 3 || s == 5 || s == 22 %only ran 180 trials
    matrixHeight = 600;
    R3 = R3(1:180);
end

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
sortedAll3 = sortrows(sortedAll,1); %based on column 1 only


%Summary for "sortedAll":
%Col 1 - image index
%Col 2 - manipulation type
%Col 3 - correct/incorrect

%For every nth row, extract reaction time and put in column 3
%key_resp8 is regular image rxn time; key_resp9 is masked image rxn time
%key_resp_8.keys_raw, or column AE is keypress
%reaction time is key_resp_8.started_raw or column AF

%key_resp_9.keys_raw, or column AM is keypress for MASKED
%reaction time is key_resp_8.started_raw or column AN for MASKED
%But hardly any are recorded.

%Maybe, first only look at Mooney recognition reaction times

%% Step 2: Add eyetracking (1 or 0) data to invalidate some trials

if subjNum ~= 22
[invalidTrial,invalidVF,VFtrials] = eyetrack_analysis3(s,edffile,T,T2,sortedManipList,sortedImgsList);
end

close all

%for subject 1 & 17: eyetracking data did not get saved properly.
if subjNum == 17 || subjNum == 22
    VFtrials = getVFtrials(sortedManipList,sortedImgsList);
    invalidVF = NaN;
end


allinvalid(s) = length(invalidVF);

if eyetrackcorrect == 'n'
    invalidVF = VFtrials;
else
    for trialnum = 1:length(sortedAll2)
        if ismember(sortedAll2(trialnum,1),invalidVF)
            sortedAll2(trialnum,3) = NaN;
        end
    end
end

%% Step 3: Analyze % correct for pre/GS/post, sorted by manipulation type

% Add correctness in column 2
totalImages = 219; %not 219 images, but image index labels range from 1-219

mooney_yesno = repmat([0; 1; 0],matrixHeight/10,1); %0 if it's a mooney, 1 if it's a GS
mooney_prepost = repmat([1;0;0],matrixHeight/10,1); %1 if it's a pre-M, 0 if it's a post-M


for gM = 1:10
    for rr = 1:totalImages
        correct.manip{gM,1}(:,subjNum) = nan(totalImages,1); %pre
        correct.manip{gM,2}(:,subjNum) = nan(totalImages,1); %grayscale
        correct.manip{gM,3}(:,subjNum) = nan(totalImages,1); %post
    end
end

imghistory = [];
for imgIdx = 1:length(sortedAll3)
    if ~isnan(sortedAll3(imgIdx,1))
        if sortedManipList(imgIdx) == 0 && ~isnan(sortedAll3(imgIdx))%if it's a Mooney image
            img = sortedImgsList(imgIdx);
            if ~ismember(img,imghistory) % if it's a Pre Mooney
                maniprows = sortedAll2(:,2);
                resprows = sortedAll2(:,3);
                gM = sum(maniprows(sortedAll2(:,1) == img));
                resp = resprows(sortedAll2(:,1) == img);
                correct.manip{gM,1}(img,s) = resp(1); % first one is pre-M resp
                correct.manip{gM,2}(img,s) = resp(2);
                correct.manip{gM,3}(img,s) = resp(3);
                imghistory = [imghistory; img];
            end
        end
    end
end


% Important: save correct.mat
save('correct_3dvaVF.mat','correct');


    end
end

%% Average & Plot:
load('C:\Users\ahach\OneDrive\Desktop\MOONEY_DATA\data_all\correct_3dvaVF.mat');

VFmanips = [1 4 10]; % Only recalculate % correct for orig + VF trials + catch
for gg = 1:3
    gM = VFmanips(gg);
for ss = 1:size(correct.manip{4,1},2)
    if ~isnan(nanmean(correct.manip{4,2}(:,ss)))
        VFsubjnum(ss) = 1;
        correctVFtrial{gg,1}(:,ss) = correct.manip{gM,1}(:,ss);
        correctVFtrial{gg,2}(:,ss) = correct.manip{gM,2}(:,ss);
        correctVFtrial{gg,3}(:,ss) = correct.manip{gM,3}(:,ss);
    else
        VFsubjnum(ss) = 0;
        correctVFtrial{gg,1}(:,ss) = nan(219,1);
        correctVFtrial{gg,2}(:,ss) =  nan(219,1);
        correctVFtrial{gg,3}(:,ss) =  nan(219,1);
    end
end
end

%mean by subject

for gg = 1:3
    for subj = 1:33 %for each subject
        for trial_type = 1:3 %3 because preM/GS/postM
            meanbysubjVF{gg,trial_type}(:,subj) = nanmean(correctVFtrial{gg,trial_type}(:,subj));
        end
    end
end

%all subject mean
for gg = 1:3
for trial_type = 1:3 %3 because preM/GS/postM
    subjAllmeanVF{gg,trial_type} = nanmean(meanbysubjVF{gg,trial_type});
end
end

%% For VF stats with 3dva inclusion radius
manipTypes = {'Original','VFshift','Catch'};
origColor = [99/255 169/255 95/255];
catchColor = [250/255, 125/255, 0];
mooneyLabel = [{'Pre'}; repmat({'Pre'}, 29, 1); {'Post'}; repmat({'Post'}, 29, 1)];
clearvars stats
stats(:,1) = meanbysubjVF{1,1}';
stats(:,2) = meanbysubjVF{1,3}';
stats(:,3) = meanbysubjVF{2,1}';
stats(:,4) = meanbysubjVF{2,3}';
stats(:,5) = meanbysubjVF{3,1}';
stats(:,6) = meanbysubjVF{3,3}';

% stats = statsall{4};
clearvars vioPts
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


