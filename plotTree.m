function plotTree(T)

adjMatrix = zeros(max(T.clu_id));

for index = 2:length(T.pa)   
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

nodeHandlers = BGH.Nodes;

 for indexNode = 1:length(nodeHandlers)
     set(BGH.Nodes(indexNode),'Size', 1+floor([plotClusterSizes(indexNode)/sizeScale plotClusterSizes(indexNode)/sizeScale]), 'Shape','circle')
 end

% g = biograph.bggui(BGH);
% %             adjustBGPlot
% f = figure('Color', 'w','pos',[520    80   554   718]);
% copyobj(g.biograph.hgAxes,f);
% delete(BGH)
% delete(g)

end