function results = reductionSimulations(results, projections)
% Date 5/14/14: written by Eugenio Marco
% Calculate perturbation effects along bifurcation axis

makePlot = 0;

%%%%Do we have to use such as high power?
reductionPowers = [1 log2(10) log2(100)]; % This means reductionFactor = 2^reductionPower;


results.projections = projections;
%results.projections = results.v';

results = getSplittingProbabilities(results, makePlot);

results.splitProbUnstableState = getInterpolatedSplitProb(results,results.roots(2));

nGenes = length(results.genes);

mutantSplits = zeros(nGenes,length(reductionPowers));

for index = 1:nGenes
    for index2 = 1:length(reductionPowers)
        reductionFactor = -results.projections(index);
        
        scoreReduced = results.roots(2) +...
            reductionPowers(index2)*reductionFactor;
        %       Sometimes we have to flip the sign
        splitProbReduced = getInterpolatedSplitProb(results,scoreReduced);
        mutantSplits(index,index2) = mean(splitProbReduced);
    end
end

results.mutantSplits = mutantSplits;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function splitProb = getInterpolatedSplitProb(results,x)

xSplit = results.xSplit;
splitProbRight = results.splitProbRight;

splitProb = zeros(size(x));

for index=1:length(x)
    if x(index)<xSplit(1)
        splitProb(index) = splitProbRight(1);
    elseif x(index)>xSplit(end)
        splitProb(index) = splitProbRight(end);
    else
        splitProb(index) = spline(xSplit,splitProbRight,x(index));
    end
end

end

function results = getSplittingProbabilities(results, makePlot)
% Calculation of Splitting probabilities

x = results.x;

unstableRoot = results.roots(2);

% expPotential = calculateExpPotential(x,shift,mu_sigma2,a,k_sigma2);
% expPotential = calculateExpPotential(x,p(1),p(2),p(3),p(4));


leftWellPosition = results.roots(1);
rightWellPosition = results.roots(3);

leftWellIndex = find(x<leftWellPosition, 1, 'last' );
rightWellIndex = find(x<rightWellPosition, 1, 'last' );

expPotential = exp(results.potential_b(leftWellIndex:rightWellIndex));

normalizationSplitProb = sum(expPotential)*(x(2)-x(1));
splitProbRight = cumsum(expPotential)*(x(2)-x(1))/...
    normalizationSplitProb;

results.leftWellIndex =leftWellIndex;
results.rightWellIndex =rightWellIndex;
results.xSplit = x(leftWellIndex:rightWellIndex);
results.splitProbRight = splitProbRight;

% unstableRootIndex = find(x<unstableRoot, 1, 'last' );

if makePlot
    subplot(223)
    xSplit = x(leftWellIndex:rightWellIndex);
    [AX,H1,H2] = plotyy(xSplit,splitProbRight,x,results.potential_b);
    set(AX(1),'Xlim',[results.plotLimits])
    set(AX(2),'Xlim',[results.plotLimits],'ylim', [-5 20],'ytick',[-5 0 5 10 15 20],...
        'yticklabel',[-5 0 5 10 15 20],'xticklabel', '')
    set(get(AX(1),'Ylabel'),'String','Splitting probability')
    set(get(AX(2),'Ylabel'),'String','Potential')
    set(H1,'LineStyle','--')
    set(H2,'LineStyle',':')
    hold on
    title('Splitting probability: falling on the right well as a function of initial state')
    line([unstableRoot unstableRoot],[0 1],'color','r')
    HHH= gcf;
end

end