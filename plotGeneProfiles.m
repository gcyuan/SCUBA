function plotGeneProfiles(dataset, selectedGenes)
% function to plot profiles of selected genes ordered by pseudotime

[~, processDataMat, ~, ~, ~, ~, ~, ~] = initialization(dataset);

load(processDataMat)

[~, locs] = ismember(selectedGenes, pro.gname);

pro.pseudotime = (pro.pseudotime- min(pro.pseudotime))/(max(pro.pseudotime)-min(pro.pseudotime));

[pseudotimeSorted, pseudotimeSortIndices] = sort(pro.pseudotime);

expr = pro.expr(pseudotimeSortIndices,locs);
   
nWindows = 100;
t0 = pseudotimeSorted(1);
t1 = pseudotimeSorted(end);

stepSize = (t1-t0)/2/nWindows;

boundariesAndCenters = t0:stepSize:t1;
windowCenters = boundariesAndCenters(2:2:end);

windowRadius = 0.08*(t1-t0);

smoothExpr = calculateWindowMedians(windowCenters, windowRadius, pseudotimeSorted, expr);


%%
figure
colorList = {'b', 'r', 'g', 'k', 'm'};
for index = 1:size(smoothExpr,2)
    plot(windowCenters, smoothExpr(:,index),colorList{index},'linewidth',2)
    hold on
end
set(gca,'fontsize',18)
legend(selectedGenes, 'location','SouthEast', 'fontsize',12)
xlabel('Pseudotime')
ylabel('Normalized Intensity')
ylim([0 1.05])
xlim([0 max(pro.pseudotime)])
title(dataset)


%%
end

function smoothExpr = calculateWindowMedians(windowCenters, windowWidths, pseudotime, expr)

smoothExpr = zeros(length(windowCenters),size(expr,2));

for index = 1:length(windowCenters)
    windowLocations = pseudotime > windowCenters(index) - windowWidths & ...
         pseudotime < windowCenters(index) + windowWidths;
     smoothExpr(index,:) = median(expr(windowLocations,:));
end

smoothExpr = smoothExpr./repmat(max(smoothExpr),length(windowCenters),1);


end
