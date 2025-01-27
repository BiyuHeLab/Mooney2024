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
