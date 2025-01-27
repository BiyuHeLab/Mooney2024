%% Computing overlapping RFs
% Neighboring voxels with poor pRF model fits (variance
% explained <50%) were removed from this computation.

clear;
clc;
close all
% cd datapath

%% Upload feature images:
% cd 'datapath\features';
percentOverlap = 0.2;
% P = 'datapath\images\features';
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
    imgFt{k,3} = imcrop(imresize(img_orig,2),[126,126,499,499]); %big
    imgFt{k,4} = [img_orig(:,251:500) zeros_matrix]; %left shift
    imgFt{k,5} = [zeros_matrix img_orig(:,1:250)];%right shift
    imgFt{k,6} = imrotate(img_orig,90); %CW rotate
    imgFt{k,7} = imrotate(img_orig,270); %CCW rotate
    imgFt{k,8} = flip(img_orig ,2); %LR invert; 2 is for horizontal flip
end

%% Load pRF values from Kay, et. al (2013):
load('grandsummary.mat');

numROI = 7;
cutoff = 50;
numSession = 9;

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

%% Load pRFs from Kay, et. al (2015):
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

%% sanity check plots

figure;
set(gcf,'color','w');
roinames = {'V1','V2','V3','V3AB','hV4','LO1','LO2','IOG','pFus','mFus'};
for roi = 1:7
subplot(3,4,roi);
abovecutoff = find(prfdata{roi}.R2>cutoff);
ecc = sqrt((prfdata{roi}.x(abovecutoff)).^2 + (prfdata{roi}.y(abovecutoff)).^2);
scatter(ecc,prfdata{roi}.sz(abovecutoff),3,'filled');
xlabel('pRF eccentricity (deg)');
ylabel('pRF size (deg)');
ylim([0 20]);
xlim([0 20]);
title(roinames{roi});
end

%% Plotting sanity check that the RF look okay.
hardcut = 50;
stimdeg = 12;

data=prfdata; %if plotting 2013 & 2015 pRFs
% data=alldata3; %if plotting 2015 pRFs only

for region = 1:length(data)

centerX = data{region}.x;
centerY = data{region}.y;
sz = data{region}.sz;

ii = data{region}.R2 >= hardcut; %R2 (model fit) cut-off
sz2 = sz(ii);
sz3 = sz2(sz2~=Inf);
sz4 = sz3(sz3<stimdeg);

centerX2 = centerX(ii);
centerX3 = centerX2(sz2~=Inf);
centerX4 = centerX3(sz3<stimdeg);

centerY2 = centerY(ii);
centerY3 = centerY2(sz2~=Inf);
centerY4 = centerY3(sz3<stimdeg);

s = centerX4 < stimdeg & centerY4 < stimdeg & centerX4 > -stimdeg & centerY4 > -stimdeg; %stimdeg cut-off to constrain to image size
centerX5 = centerX4(s);
centerY5 = centerY4(s);
sz5 = sz4(s);

totalRF = length(sz5);

figure(region)
hold on
    for i = 1:20%totalRF
        rr = sz5(i)*2;
        xy = [centerX5(i);centerY5(i)];
    % https://www.mathworks.com/matlabcentral/answers/67605-draw-circle-on-image

    %         figure;
        c = xy';
        pos = [c-rr 2*rr 2*rr];
        r = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'None', 'Edgecolor','w');
        axis equal
        axis off
        title('All RF plot, sanity check (in DVA)');
    end
end

%% Now, quantify overlap between RF + feature
featureRF = cell(length(prfdata),length(imgFt));
hardcut = 50;
stimdeg = 24;

% maxStimDeg = [1.8 2.1 2.2 4.1 8.1 8.2 13];
for region = 1:length(prfdata)
    disp('Region'); disp(region);
   
    centerX = prfdata{region}.x;
    centerY = prfdata{region}.y;
    sz = prfdata{region}.sz;

    ii = prfdata{region}.R2 >= hardcut; %R2 (model fit) cut-off
    sz2 = sz(ii);
    sz3 = sz2(sz2~=Inf);
    sz4 = sz3(sz3<stimdeg);

    centerX2 = centerX(ii);
    centerX3 = centerX2(sz2~=Inf);
    centerX4 = centerX3(sz3<stimdeg);

    centerY2 = centerY(ii);
    centerY3 = centerY2(sz2~=Inf);
    centerY4 = centerY3(sz3<stimdeg);

    s = centerX4 < stimdeg & centerY4 < stimdeg & centerX4 > -stimdeg & centerY4 > -stimdeg; %stimdeg cut-off to constrain to image size
    centerX5 = centerX4(s);
    centerY5 = centerY4(s);
    sz5 = sz4(s);

    totalRF = length(sz5);

    for i = 1:totalRF

        rr = sz5(i)*2;
        xy = [centerX5(i);centerY5(i)];
        imageSizeX = 500; %These meshgrid pixel boundaries exclude RF outside image
        imageSizeY = 500;
        [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

    % https://www.mathworks.com/matlabcentral/answers/67605-draw-circle-on-image
        c = xy';
        centerX = c(1)*(500/12)+250;
        centerY = c(2)*(500/12)+250;
        radius = rr*(500/12)/2;
        circlePixels = (rowsInImage - centerY).^2 ...
            + (columnsInImage - centerX).^2 <= radius.^2;

    for jj = 1:length(imgFt)
            %if the circle overlaps with feature, mark as 1
            RFoverlap = imgFt{jj,1} & circlePixels;
            if sum(sum(RFoverlap)) > 1
                featureRF{region,jj}(i,1) = 1;
                for m = 1:7 %if there is an overlap with this RF, check if there is an overlap bw RF and all img manipulated features.
                RFoverlap2 = imgFt{jj,m+1} & circlePixels;
                    if sum(sum(RFoverlap2))/sum(sum(RFoverlap)) > percentOverlap
                        featureRF{region,jj}(i,m+1) = 1;
                    else
                        featureRF{region,jj}(i,m+1) = 0;
                    end
                end
           
            else
                featureRF{region,jj}(i,1) = 0;
                featureRF{region,jj}(i,2:8) = 0;
            end
    end
    end

end

%% Count the number of prfdata RFs.
for region = 1:length(prfdata)
    for jj = 1:length(imgFt)
        for m = 1:8
            tempRF(jj,m) = length(find(featureRF{region,jj}(:,m) == 1))/length(find(featureRF{region,jj}(:,1)==1));
        end
    end
        overlappingRF(:,region) = mean(tempRF);
end

color = {'k',[0.7 0 1],[0.8 0.6 0.8],'b',[0.3010 0.7450 0.9330],'g',[0.4660 0.6740 0.1880],'r'};

%delete V3AB:
overlappingRF(:,4) = [];

figure; hold on
set(gcf,'color','w');
for m = 1:8
    p(m) = plot(overlappingRF(m,:)*100,'Color',color{m},'LineWidth',2);
    title('% Overlapping RF across visual hierarchy');
    xlim([0 11]);
    ylim([-10, 110]);
    fname = ({'V1','V2','V3','hV4','LO1','LO2','IOG','pFus','mFus'});
    set(gca, 'XTick', 1:length(fname),'XTickLabel',fname,'FontSize',18);
    ylabel('% Overlapping RFs');
    xlabel('ROIs')
end
legend(manipTypes,'Location','southeastoutside');
ax = gca;
ax.FontSize = 11;
legend(manipTypes);

