clear;
clc;
close all
cd 'C:\Users\ahach\OneDrive\Desktop\Mooney Project\InPersonFiles\images\features'
%% Upload feature images:
% cd 'C:\Users\ahach\OneDrive\Desktop\InPersonFiles\images\features';
percentOverlap = 0.25;

P = 'C:\Users\ahach\OneDrive\Desktop\Mooney Project\InPersonFiles\images\features';
S = dir(fullfile(P,'*.bmp'));
imgFt = cell(size(S,1),1);
for k = 1:numel(S)
    F = S(k).name; 
    temp = imread(F); %load image
    imgFt{k} = logical(double(temp)); %convert to double
end

%% Manipulate feature images:
zeros_matrix = zeros(500,250);
manipTypes = {'Original','Small','Big','Lshift','Rshift','CWrot','CCWrot','LRinvert'};

for k = 1:numel(S)
    img_orig = imgFt{k,1};
    imgFt{k,2} = padarray(imresize(img_orig,0.5),[125 125],0,'both'); %small
    imgFt{k,3} = imresize(img_orig,2,'nearest'); %big
    imgFt{k,4} = [img_orig(:,251:500) zeros_matrix]; %left shift [img_orig zeros_matrix]; %left shift
    imgFt{k,5} = [zeros_matrix img_orig(:,1:250)];%right shift [zeros_matrix img_orig];%right shift
    imgFt{k,6} = imrotate(img_orig,90); %CW rotate
    imgFt{k,7} = imrotate(img_orig,270); %CCW rotate
    imgFt{k,8} = flip(img_orig ,2); %LR invert; 2 is for horizontal flip
end
% for k = 1:numel(S)
%     img_orig = imgFt{k,1};
%     imgFt{k,2} = padarray(imresize(img_orig,0.5),[125 125],0,'both'); %small
%     imgFt{k,3} = imcrop(imresize(img_orig,2),[126,126,499,499]); %big
%     imgFt{k,4} = [img_orig(:,251:500) zeros_matrix]; %left shift
%     imgFt{k,5} = [zeros_matrix img_orig(:,1:250)];%right shift
%     imgFt{k,6} = imrotate(img_orig,90); %CW rotate
%     imgFt{k,7} = imrotate(img_orig,270); %CCW rotate
%     imgFt{k,8} = flip(img_orig ,2); %LR invert; 2 is for horizontal flip
% end

% for k = 1:numel(S)
%     figure;
%     imshow(imgFt{k,3});
% end

%% Load pRF values from Kendrick (2013):
load('grandsummary.mat');

numROI = 7;
cutoff = 50;
numSession = 9;

%For all sessions:
% prfdata = cell(numROI,1);
% for s = 1:numSession
%     for r = 1:numROI
%         roiInd = find(wtf(s).roi == r);
%         for rr = 1:length(roiInd)
%             roi = roiInd(rr);
%             prfdata{r}.x(rr) = wtf(s).x(roi);
%             prfdata{r}.y(rr) = wtf(s).y(roi);
%             prfdata{r}.sz(rr) = wtf(s).sigma(roi);
%             prfdata{r}.R2(rr) = wtf(s).rnonlinear(roi);
%         end
%     end
% end

%For contrast patterns:
prfdata = cell(numROI,1);
sess = [3 5 7];
earlyROI = 1:5;
for i = 1:length(sess)
    s = sess(i);
    for j = 1:length(earlyROI)
        r = earlyROI(j);
        roiInd = find(wtf(s).roi == r & wtf(s).isin == 1);
        for rr = 1:length(roiInd)
            roi = roiInd(rr);
            prfdata{j}.x(rr) = wtf(s).x(roi);
            prfdata{j}.y(rr) = wtf(s).y(roi);
            prfdata{j}.sz(rr) = wtf(s).sigma(roi);
            prfdata{j}.R2(rr) = wtf(s).rnonlinear(roi);
        end
    end
end

%For object patterns:
sess = [1 2 4 5 8 9];
objROI = [8 9];
for i = 1:length(sess)
    s = sess(i);
    for j = 1:length(objROI)
        r = objROI(j);
        roiInd = find(wtf(s).roi == r & wtf(s).isin == 1);
        for rr = 1:length(roiInd)
            roi = roiInd(rr);
            prfdata{j+length(earlyROI)}.x(rr) = wtf(s).x(roi);
            prfdata{j+length(earlyROI)}.y(rr) = wtf(s).y(roi);
            prfdata{j+length(earlyROI)}.sz(rr) = wtf(s).sigma(roi);
            prfdata{j+length(earlyROI)}.R2(rr) = wtf(s).rnonlinear(roi);
        end
    end
end

%% Load pRFs from 2015 paper:
load('kaycurbio2.mat');
allRegions = length(alldata);

% https://www.mathworks.com/matlabcentral/answers/96973-how-can-i-concatenate-or-merge-two-structures
alldata2{1,1} = MergeStructs(alldata{1,1},alldata{1,2});
alldata2{1,2} = MergeStructs(alldata{1,3},alldata{1,4});
alldata2{1,3} = MergeStructs(alldata{1,5},alldata{1,6});
alldata2{1,4} = MergeStructs(alldata{1,7},alldata{1,8});
alldata2{1,5} = MergeStructs(alldata{1,9},alldata{1,10});
alldata2{1,6} = MergeStructs(alldata{1,11},alldata{1,12});
alldata2{1,7} = MergeStructs(alldata{1,13},alldata{1,14});

for i = 1:length(alldata2)
    temp = alldata2{1,i};
    alldata3{1,i}.R2 = [temp(1).R2; temp(2).R2];
    alldata3{1,i}.x = [temp(1).x; temp(2).x];
    alldata3{1,i}.y = [temp(1).y; temp(2).y];
    alldata3{1,i}.sz = [temp(1).sz; temp(2).sz];
end

prfdata{8} = alldata3{5}; %IOG
prfdata{9} = alldata3{6}; %pFus
prfdata{10} = alldata3{7}; %mFus

%% Plotting pRFs & images:

%counting pRFs per ROI:
hardcut = 50;
for r = 1:length(prfdata)
    pRFcount(r) = length(prfdata{r}.R2(prfdata{r}.R2 >= hardcut)); %R2 (model fit) cut-off
end

hardcut = 80;
stimdeg = 12/4;
ones_matrix = ones(500,250)*255;
data=prfdata; %if plotting 2013 pRFs
% data=alldata3; %if plotting 2015 pRFs
img = imread('C:\Users\ahach\OneDrive\Desktop\Manuscript\1048.bmp');
imgsmall = padarray(imresize(img,0.5),[125 125],255,'both'); %small
imgbig = imresize(img,2,'nearest'); %big
imgleft = [img ones_matrix]; %left shift
imgright = [ones_matrix img];%right shift
imgCW = imrotate(img,90); %CW rotate
imgCCW = imrotate(img,270); %CCW rotate
imgLR = flip(img ,2); %LR invert; 2 is for horizontal flip

[rows, cols] = size(img);
figure;
counter=1;
for region = [1 2 3 5 6 7 8 9 10]%1:length(data)
    for m = 1:7
        img = imread('C:\Users\ahach\OneDrive\Desktop\Manuscript\1048.bmp');
        centx = rows/2;
        centy = cols/2;
        if m == 1
            img = img;
            centx = rows/2;
            centy = cols/2;
        elseif m == 2
            img = imgLR;
            centx = rows/2;
            centy = cols/2;
        elseif m == 3
            img = imgCW;
            centx = rows/2;
            centy = cols/2;
        elseif m == 4
            img = imgCCW;
            centx = rows/2;
            centy = cols/2;
        elseif m == 5
            img = imgleft;
            centx = rows/2 + 250;
            centy = cols/2;
        elseif m == 6
            img = imgright;
            centx = rows/2;
            centy = cols/2;
        elseif m == 7
            img = imgsmall;
            centx = rows/2;
            centy = cols/2;
        end
    
    % aa = subplot(10,7,counter);
    subaxis(9,7,counter,'Spacing',0.01,'Margin',0);
    
    hold on
    
    imshow(img);
    
    centerX = data{region}.x;
    centerY = data{region}.y;
    sz = data{region}.sz;
    
    ii = (data{region}.R2 >= hardcut); %R2 (model fit) cut-off
    sz2 = sz(ii);
    sz3 = sz2(sz2~=Inf);
    sz4 = sz3(sz3<stimdeg);
    
    centerX2 = centerX(ii);
    centerX3 = centerX2(sz2~=Inf);
    centerX4 = centerX3(centerX3<stimdeg/2 & centerX3 > -stimdeg/2);
    
    centerY2 = centerY(ii);
    centerY3 = centerY2(sz2~=Inf);
    centerY4 = centerY3(centerY3<stimdeg/2 & centerY3 > -stimdeg/2);
    
    % s = centerX4 < stimdeg & centerY4 < stimdeg & centerX4 > -stimdeg & centerY4 > -stimdeg; %stimdeg cut-off to constrain to image size
    centerX5 = centerX4;
    centerY5 = centerY4;
    % centerX5 = centerX3;
    % centerY5 = centerY3;
    sz5 = sz4;
    
    totalRF = length(sz5);
    
        for i = 1:15%totalRF
            sz5 = rmoutliers(sz5);
            rr = sz5(i)*2*(500/12); %sz5 is in radius, so multiply by 2 to get diameter. Then, multiply by (500/12) to get units into pixels.
            xy = [centerX5(i)*(500/12) + centx;centerY5(i)*(500/12) + centy];
        % https://www.mathworks.com/matlabcentral/answers/67605-draw-circle-on-image
            c = xy';
            pos = [c-rr 2*rr 2*rr];
            r = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'None', 'Edgecolor',[255/255 189/255 31/255]);
            axis equal
            axis off
            %title('All RF plot, sanity check (in DVA)');
        end
        counter=counter+1;
    end
end
