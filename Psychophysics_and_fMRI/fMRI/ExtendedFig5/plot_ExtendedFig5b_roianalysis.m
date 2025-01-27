figure;set(gcf,'color','w');
numROIs = 19;
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

