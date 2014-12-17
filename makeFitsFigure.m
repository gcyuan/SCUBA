function makeFitsFigure(results, dataFolder)
% Date 5/14/14: written by Eugenio Marco
% Display fits

nBins = 20;

[nFinal xout] = hist([results.scoreFinalCluster1; results.scoreFinalCluster2], nBins);
nFinal1 = hist(results.scoreFinalCluster1,xout);
nFinal2 = hist(results.scoreFinalCluster2,xout);
nParent = hist(results.scoreParentCluster,xout);
normFinal = sum(nFinal)*(xout(2)-xout(1));
p = results.p;
x = results.x;
limits = results.limits;
ymleFinal = normFinal*here_steadyStateDistribution(x,p(1),p(2),p(end), limits);
normParent = sum(nParent)*(xout(2)-xout(1));
ymleParent = normParent*here_steadyStateDistribution(x,p(1),p(2),p(end-1), limits);

switch results.options.parentLevelFit
    case 'parent'
        nrows = 2;
        potential(1) = {results.potential_p};
        potential(2) = {results.potential_b};
        nData(1) = {nParent};
        nData(2) = {nFinal1};
        nData(3) = {nFinal2};
        mlePrediction(1) = {ymleParent};
        mlePrediction(2) = {ymleFinal};
    case 'parentParent'
        nrows = 3;
        potential(1) = {results.potential_pp};
        potential(2) = {results.potential_p};
        potential(3) = {results.potential_b};
        nParentParent = hist(results.scoreParentParentCluster,xout);
        normParentParent = sum(nParentParent)*(xout(2)-xout(1));
        ymleParentParent = normParentParent*here_steadyStateDistribution(x,p(1),p(2),p(3), limits);
        nData(1) = {nParentParent};
        nData(2) = {nParent};
        nData(3) = {nFinal1};
        nData(4) = {nFinal2};
        mlePrediction(1) = {ymleParentParent};
        mlePrediction(2) = {ymleParent};
        mlePrediction(3) = {ymleFinal};
end
%%
ylimits = [(min(cell2mat(potential))-3) 10];
figure
nplot = 1;
for index = 1:nrows

    subplot(nrows,2,nplot)
    plot(results.x,potential{index})
    nplot = nplot+1;
    xlim(results.plotLimits)
    title('Potential')
    xlabel('Scores on Bifurcation Axes')
    ylim(ylimits)
    title(['\sigma = ' num2str(p(1)) '  b = ' num2str(p(2))...
        '  a = ' num2str(p(index+2)) ])    
    subplot(nrows,2,nplot)
    if index<nrows
        H = bar(xout,nData{index});
    else
        H = bar(xout,[nData{index}; nData{index+1}]',2);
    end
    hold on
    plot(x,mlePrediction{index},'r')
    if index<nrows
        set(H(1),'facecolor','blue') % use color name
    else
        set(H(1),'facecolor','red') % or use RGB triple
        set(H(2),'facecolor','blue') % use color name
    end
    xlabel('Scores on Bifurcation Axes')
    ylabel('Counts')
    xlim(results.plotLimits)
    title('PDF')
    nplot = nplot+1;
end

HA = annotation('textbox',[.4 .025 .35 .05]);

set(HA,'string', ['Stage ', int2str(results.stage), ' bifurcation '],'Fontsize',18,'EdgeColor','none')


filenameFigure = fullfile(dataFolder,'figures',['results_fit_' ...
    num2str(results.stage), '.jpg']);
saveas(gcf, filenameFigure,'jpg')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ssDist = here_steadyStateDistribution(x,sigma,b,a, limits)

lowerLim = -limits;
upperLim = limits;

% Normalization
myFun = @(y) exp((y.^2.*(a-0.5*y.^2)+2*b*y)/sigma^2);
if verLessThan('matlab','8.0.1') % older versions do not have integral
    normalization = 1/quad(myFun,lowerLim,upperLim);
else
    normalization = 1/integral(myFun,lowerLim,upperLim);
end
rest = myFun(x);

ssDist = normalization*rest+realmin;


