function plotTree(T, dataset)

adjMatrix = zeros(max(T.clu_id));

nZeros = sum(T.pa==0);

for index = nZeros+1:length(T.pa)   
    adjMatrix(T.pa(index),T.clu_id(index)) = 1;
end


clusterSizes = zeros(max(T.s),1);

for index = 1:length(clusterSizes)
    clusterSizes(index) = sum(T.s==index);    
end

plotClusterSizes = clusterSizes/max(clusterSizes)*50;

sizeScale = 3;

a = 1:max(T.clu_id);
IDS = cellstr(num2str(a(:)));

BG = biograph(adjMatrix, IDS);
BGH = view(BG);

% %From the old version
% nodeHandlers = BGH.Nodes;
% 
%  for indexNode = 1:length(nodeHandlers)
%      set(BGH.Nodes(indexNode),'Size', 1+floor([plotClusterSizes(indexNode)/sizeScale plotClusterSizes(indexNode)/sizeScale]), 'Shape','circle')
%  end
% 
% g = biograph.bggui(BGH);
% %             adjustBGPlot
% f = figure('Color', 'w','pos',[520    80   554   718]);

% %Replace by the newer version 1
% %             adjustBGPlot
% f = figure('Color', 'w', 'pos',g.biograph.hgAxes.Position .* [1 1 1.5 1.5]);
% copyobj(g.biograph.hgAxes,f);

% %Replace by the newer version 2
g = biograph.bggui(BGH);
f = figure('Color', 'w', 'pos',g.biograph.hgAxes.Position .* [1 1 1.5 1.5]);
copyobj(g.biograph.hgAxes,f);

delete(BGH)
delete(g)

[~, ~, ~, ~, ~, ~, ~, figuresDir] = initialization(dataset);


filenameFigure = fullfile(figuresDir, 'tree.jpg');
saveas(gcf, filenameFigure,'jpg')

end