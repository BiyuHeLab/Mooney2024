function [invalidTrial, invalidVF,VFtrials] = eyetrack_analysis3(s,edffile,T,T2,sortedManipList,sortedImgsList)

addpath('toolboxes/edf-converter-master/@Edf2Mat/');
edf1 = Edf2Mat(edffile);
VFtrials = getVFtrials(sortedManipList,sortedImgsList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixel_per_deg = 43.6;
center_x = 960;
center_y = 540;
r = pixel_per_deg*3; %110.93 pixels is 2 dva

edfTime = edf1.Events.Messages.time;
edfEvent = edf1.Events.Messages.info;
edfTimeAll = edf1.Samples.time; %not just the time that corresponds to events var, but time for each time stamp

temp = zeros(length(edfEvent),1);

for i = 1:length(edfEvent)
    if contains(edfEvent{i},'stim_onset')
        temp (i)= edfTime(i);
    end
end

stimTime = temp(temp~=0);

stimTime_startend(:,1) = stimTime(1:4:end);
stimTime_startend(:,2) = stimTime(4:4:end);
%annoyingly, every 4, because eyetracking records masking trials that never
%happened.

%col 1 is stimulus onset, col 2 is stimulus offset

%From plot.mat
if ~exist('startIdx', 'var')
    startIdx = 1;
end

if ~exist('endIdx', 'var')
    endIdx = size(edf1.Samples.posX, 1);
end

range = startIdx:endIdx;

posX = edf1.Samples.posX(range);
% Y must be inverted, because eyetracker origin
% is upper left corner in a graph its the lower left
posY = edf1.Samples.posY(range);

%Take position x,y from stimulus 1 onset to stimulus 1 offset
gazeXbyTrial = cell(length(VFtrials),1);
gazeYbyTrial = cell(length(VFtrials),1);
for trial = 1:length(VFtrials)
    trialInd = VFtrials(trial);
    stimStart = stimTime_startend(trialInd,1); %/1000; %divide by 1000 because recording is 1000Hz; posX and posY have ~1268581
    %values. The last stimulation onset time is 9847191
    %actually no need to divide by 1000, just know that pos and time
    %sampling rate is the same.
    stimEnd = stimTime_startend(trialInd,2); %/1000;
    
    stimStartIdx = find(edfTimeAll == stimStart);
    stimEndIdx = find(edfTimeAll == stimEnd);
    gazeXbyTrial {trial} = posX(stimStartIdx:stimEndIdx);
    gazeYbyTrial {trial} = posY(stimStartIdx:stimEndIdx);
    
    subplot(3,3,trial)
    hold on
    plot(posX(stimStartIdx:stimEndIdx)/pixel_per_deg, posY(stimStartIdx:stimEndIdx)/pixel_per_deg,  '*', 'MarkerSize',2);
    title('Plot of the eye movement');
    % axis([0 1000 0 1920]);
%     axis('square');
    xlabel('x (dva)');
    ylabel('y (dva)');
    ylim([0 24]);
    xlim([0 24]);

    %43 pixels per degree, on average
    %6 degrees is 43*6 = 258
    %To be conservative, let's say 2 degree radius circle (RF size of V4).
    %That means, even the farthest right or left that they look, radius is
    %only 4 degrees from fixation.
    %43*2 = 86
    %draw a 6x6 degree circle:
    h = draw_circle(center_x/pixel_per_deg,center_y/pixel_per_deg,r/pixel_per_deg);
    axis equal
end

saveas(gcf, fullfile([pwd '/eyetrackPlots_3dvaVF/'], [edffile '.jpg']));

%% How many pixels per visual angle degree?

%Screen size is 1024 x 720
%Screen width is 30cm, distance to screen is 65cm

%Key facts:
% X Resolution,
% Y Resolution:max
% Angular resolution at current gaze position (in screen pixels
% per visual degree)

%43 pixels per visual degree 

%% Final trials to be used:
%x = 960; y = 540; r = 86;

validTrial = zeros(length(VFtrials),1);
for trial = 1:length(VFtrials)
    for dataPT = 1:size(gazeXbyTrial{trial})
        x = gazeXbyTrial{trial}(dataPT); y = gazeYbyTrial{trial}(dataPT);
        square_dist = (center_x - x)^2 + (center_y - y)^2;
            if square_dist <= r^2
                valid_temp (dataPT) = 0; %if the gaze point is inside the circle
            else
                valid_temp (dataPT) = 1; %1 is bad, violates allowable gaze shift
            end
    end
            if sum(valid_temp) == 0 %if all the gaze points are inside the circle
                invalidTrial (trial) = 0;
            else
                invalidTrial (trial) = 1;
            end
end

invalidVF = zeros(length(VFtrials),1);
for trial = 1:length(VFtrials)
    if invalidTrial(trial) == 1 %if the trial is invalid because of gaze shift
        invalidVF (trial) = VFtrials(trial);
    else
        invalidVF (trial) = 0;
    end
end

invalidVF = invalidVF(invalidVF~=0);
end
