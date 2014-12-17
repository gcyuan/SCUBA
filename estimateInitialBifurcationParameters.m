function initialParameters = estimateInitialBifurcationParameters(results)

g3 = results.scoreFinalCluster1;
g4 = results.scoreFinalCluster2;


aBifurcation = ((mean(g4)-mean(g3))/2)^2;
aParent = -aBifurcation;
aParentParent = -2*aBifurcation;
b = 0;
sigma = var([g3; g4]);

initialParameters = [sigma, b , aParentParent, aParent, aBifurcation];


end