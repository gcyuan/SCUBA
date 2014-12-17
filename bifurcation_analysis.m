function bifurcation_analysis(dataset)
% Date 5/13/14: written by Eugenio Marco
%
% This code takes sthe gene expression data for several cell stages projected
% along a bifurcation axis and fits a model based on the solution of the
% Fokker-Planck equation with a quartic potential

% option for setting the depth of the fit
options.parentLevelFit = 'parentParent';
% options.parentLevelFit = 'parent';

% option to calculate the log likelihood with the average number of cells
% per cell stage
options.normalizeLikelihoodLevelCellCounts = 1;
% options.normalizeLikelihoodLevelCellCounts = 0;

[~, ~, ~, PCAdataFile, dataFolder, ~, intermediate_filesDir, figuresDir] = initialization(dataset);

treefile = fullfile(intermediate_filesDir,'final_tree.mat');

load(PCAdataFile);

load(treefile);

minNCellsInBifurcation = 30;

bif_in = fullfile(intermediate_filesDir,'bifurcation_direction.mat');
if ~exist(bif_in, 'file')
    warning('No bifurcations found');
else
    load(bif_in);
    nbif = length(bif.t);
    
    makePlot = 0;
    
    allResults = cell(1, nbif);
    
    for index = 1:nbif,
        results = getScores(bif, T, index, makePlot);
        results.genes = pro.gname';
        if length(results.scoreFinal) > minNCellsInBifurcation
            results = fitDistribution(results,dataFolder,makePlot,options);
            results = reductionSimulations(results, bif.direction(:, index));
            allResults(index) = {results};
        else
            disp(['Not enough cells to fit bifurcation potential at stage ', num2str(bif.stage(index)), '. At least ' num2str(minNCellsInBifurcation) ' cells are required.'])
        end
    end
    
    if length(results.scoreFinal) > minNCellsInBifurcation
        saveAllResults(dataFolder, allResults)
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = getScores(bif, T, index, makePlot)

parentParentCluster = bif.clu_id{index}(1);
parentCluster = bif.clu_id{index}(2);
finalCluster1 = bif.clu_id{index}(3);
finalCluster2 = bif.clu_id{index}(4);

allScores = bif.proj(:, index);


dataParentParentCluster = allScores(T.s==parentParentCluster);
dataParentCluster = allScores(T.s==parentCluster);
dataFinalCluster1 = allScores(T.s==finalCluster1);
dataFinalCluster2 = allScores(T.s==finalCluster2);

%%%%%%%Need to think about how to define results.bifurcation

% Scores recentered
clear results
results.scoreParentParentCluster = dataParentParentCluster-mean(dataParentParentCluster);
results.scoreParentCluster = dataParentCluster- mean(dataParentCluster);

centerFinal = (mean(dataFinalCluster1)+ ...
    mean(dataFinalCluster2))/2; % midpoint between final clusters

results.scoreFinal = [dataFinalCluster1; dataFinalCluster2]-centerFinal;
results.scoreFinalCluster1 = dataFinalCluster1-centerFinal;
results.scoreFinalCluster2 = dataFinalCluster2-centerFinal;

[~, xout] = hist([dataFinalCluster1-centerFinal; dataFinalCluster2-centerFinal], 20);

results.limits = max(abs(results.scoreFinal))*1.5;
results.plotLimits = [(xout(1)-0.2*(xout(end)-xout(1))) (xout(end)+0.2*(xout(end)-xout(1)))];

%results.v = v;

results.t = T.t(parentCluster);
results.stage = T.stage(parentCluster);

if makePlot
    figure
    
    nFinal1 = hist(dataFinalCluster1,xout);
    nFinal2 = hist(dataFinalCluster2,xout);
    nParent = hist(dataParentCluster,xout);
    nParentParent = hist(dataParentParentCluster,xout);
    H = bar(xout,[nParentParent; nParent; nFinal1; nFinal2]',2);
    set(H(1),'facecolor','green') % or use RGB triple
    set(H(2),'facecolor','magenta') % use color name
    set(H(3),'facecolor','red') % or use RGB triple
    set(H(4),'facecolor','blue') % use color name
    legend(['Cluster ' num2str(parentParentCluster)], ['Cluster ' num2str(parentCluster)],...
        ['Cluster ' num2str(finalCluster1)], ['Cluster ' num2str(finalCluster2)])
    xlabel('Scores on Bifurcation Axes')
    ylabel('Counts')
    xlim(results.plotLimits)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = fitDistribution(results,dataFolder,makePlot,options)

% We read tree if we want initial estimates of shifts
% mixtureTree = readMixtureTreeData;
% switch results.bifurcation(1)
%     case 32
%         InitialShifts = [mixtureTree.X32(results.bifurcation(2)) ...
%             mixtureTree.X32(results.bifurcation(3)) ...
%             mean(mixtureTree.X32(results.bifurcation(4:5)))];
%     case 64
%         InitialShifts = [mixtureTree.X64(results.bifurcation(2)) ...
%             mixtureTree.X64(results.bifurcation(3)) ...
%             mean(mixtureTree.X64(results.bifurcation(4:5)))];
% end

[~, dataset, ~] = fileparts(dataFolder);

% Sizes of parentParent, parent and after bifurcation
switch options.parentLevelFit
    case 'parentParent'
        
        data = [results.scoreParentParentCluster; results.scoreParentCluster;...
            results.scoreFinal];
        dataSizes = [length(results.scoreParentParentCluster)...
            length(results.scoreParentCluster) length(results.scoreFinal)];
        % format of initialParams is [sigma b a_pp a_p a_b]
        initialParams = estimateInitialBifurcationParameters(results);
end

x = linspace(-results.limits,results.limits,5000);

oldOptions = statset('mlecustom');
newOptions = statset(oldOptions,'MaxIter',10000,'MaxFunEvals',10000);


% disp(['Begin of the MLE to find the parameters of the ' num2str(results.bifurcation(1))...
%     '-cell bifurcation'])
% Our default
p = mle(data,'pdf',{@multiStateDistribution results.limits dataSizes options},'options',...
    newOptions,'start',initialParams);

likelihood.options = options;
likelihood.options.normalizeLikelihoodLevelCellCounts = 0;

pcell = num2cell(p);
results.likelihood = prod(multiStateDistribution(data,pcell{:},results.limits,dataSizes,likelihood.options));


results.x = x;
results.p = p;
results.initialParams = initialParams;
results.options = options;
results.dataSizes = dataSizes;

results = getPotentialDividedBySigmaSquared(results);
results = getRootsDerivativeBifurcationPotential(results);



if makePlot
    figure
    [nFinal, xout] = hist(results.scoreFinal,20);
    normFinal = sum(nFinal)*(xout(2)-xout(1));
    ymleFinal = normFinal*steadyStateDistribution(x,p(1),p(2),p(end),results.limits);
    
    nFinal1 = hist(results.scoreFinalCluster1,xout);
    nFinal2 = hist(results.scoreFinalCluster2,xout);
    nParent = hist(results.scoreParentCluster,xout);
    normParent = sum(nParent)*(xout(2)-xout(1));
    ymleParent = normParent*steadyStateDistribution(x,p(1),p(2),p(end-1),results.limits);
    H = bar(xout,[nFinal1; nFinal2; nParent]',2);
    set(H(1),'facecolor','blue') % use color name
    set(H(2),'facecolor','magenta') % use color name
    set(H(3),'facecolor','red') % or use RGB triple
    hold on
    plot(x,ymleFinal,'r')
    plot(x,ymleParent,'k')
    title(['\sigma^2 = ' num2str(p(1)) '  b = ' num2str(p(2)) '  shift = ' num2str(p(end-1)) ...
        '  a = ' num2str(p(end)) ])
    xlim(results.plotLimits)
    
    filenameFigure = fullfile(figuresDir,'fit.jpg');
    saveas(gcf, filenameFigure,'jpg')
end

makeFitsFigure(results, dataFolder)

filenameResults = fullfile(dataFolder,'intermediate_files',['results_fit_' ...
    num2str(results.stage) '.mat']);
save(filenameResults, 'results')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msDist = multiStateDistribution(data,varargin)

switch nargin
    case 9 % We go back to parentParent level
        [sigma, b, a_pp, a_p, a_b] = varargin{1:end-3};
    case 7 % We go back to parent level
        [sigma, b, a_p, a_b] = varargin{1:end-3};
end

options = varargin{end};
dataSizes = varargin{end-1};
limits = varargin{end-2};

cumsumDataSizes = cumsum(dataSizes);

splitData = cell(size(dataSizes));
splitData(1) = {data(1:dataSizes(1))};
for index = 2:length(dataSizes)
    splitData(index) = {data(cumsumDataSizes(index-1)+1:cumsumDataSizes(index))};
end

switch nargin
    case 9 % We go back to parentParent level
        dist_pp = steadyStateDistribution(splitData{1},sigma,b,a_pp,limits);
        dist_p = steadyStateDistribution(splitData{2},sigma,b,a_p,limits);
        dist_b = steadyStateDistribution(splitData{3},sigma,b,a_b,limits);
        if options.normalizeLikelihoodLevelCellCounts
            msDist = [dist_pp.^(1/dataSizes(1)); dist_p.^(1/dataSizes(2)); dist_b.^(1/dataSizes(3))];
        else
            msDist = [dist_pp; dist_p; dist_b];
        end
    case 7  % We go back to parent level
        dist_p = steadyStateDistribution(splitData{1},sigma,b,a_p,limits);
        dist_b = steadyStateDistribution(splitData{2},sigma,b,a_b,limits);
        if options.normalizeLikelihoodLevelCellCounts
            msDist = [dist_p.^(1/dataSizes(1)); dist_b.^(1/dataSizes(2))];
        else
            msDist = [dist_p; dist_b];
        end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ssDist = steadyStateDistribution(x,sigma,b,a, limits)

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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = getRootsDerivativeBifurcationPotential(results)

switch results.options.parentLevelFit
    case 'parentParent'
        a = results.p(5);
        b = results.p(2);
    case 'parent'
        a = results.p(4);
        b = results.p(2);
end

p = [-1 0 a b];

r = roots(p);

results.roots = sort(real(r));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = getPotentialDividedBySigmaSquared(results)

x = results.x;

myPotential = @(y,sigma,b,a) -(((y).^2.*(a-0.5*(y).^2) + 2*b*(y))/sigma^2);

switch results.options.parentLevelFit
    case 'parentParent'
        results.potential_pp = myPotential(x,results.p(1),results.p(2),results.p(3));
        results.potential_p = myPotential(x,results.p(1),results.p(2),results.p(4));
        results.potential_b = myPotential(x,results.p(1),results.p(2),results.p(5));
    case 'parent'
        results.potential_p = myPotential(x,results.p(1),results.p(2),results.p(3));
        results.potential_b = myPotential(x,results.p(1),results.p(2),results.p(4));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveAllResults(dataFolder, allResults)

%%
% % save reduction results in txt format
% fileReduction{1} = fullfile(dataFolder,'results','results_perturbation_32.txt');
%
% fileReduction{2} = fullfile(dataFolder,'results','results_perturbation_64.txt');
%
% textReductions = {'Reduction X32' 'Reduction X64'};
nbif = length(allResults);
filePerturbation = fullfile(dataFolder,'results','results_perturbation.txt');
fs = fopen(filePerturbation, 'w+');

headers = cell(1, nbif+1);
nheader = length(headers);
headers{1} = 'Stage';
for k = 1:nbif,
    headers{k+1} = ['Stage_', num2str(allResults{k}.stage)];
end
for k = 1:nheader-1,
    fprintf(fs, '%s\t', headers{k});
end
fprintf(fs, '%s\n', headers{end});

%%%%%% No sorting.
geneName = allResults{1}.genes;
ngene = length(geneName);
for n = 1:ngene,
    fprintf(fs, '%s\t', geneName{n});
    for index = 1:nbif-1
        genePeturbation = allResults{index}.mutantSplits-allResults{index}.splitProbUnstableState;
        fprintf(fs, '%f\t', genePeturbation(n,1));
    end
    genePeturbation = allResults{nbif}.mutantSplits-allResults{nbif}.splitProbUnstableState;
    fprintf(fs, '%f\n', genePeturbation(n,1));
end
fclose(fs);

%%

% save model parameters
fileParameters = fullfile(dataFolder,'results','fit_parameters.txt');

fs = fopen(fileParameters, 'w+');

headers = cell(1, nbif+1);
nheader = length(headers);
headers{1} = 'Stage';
for k = 1:nbif,
    headers{k+1} = ['Stage_', num2str(allResults{k}.stage)];
end
for k = 1:nheader-1,
    fprintf(fs, '%s\t', headers{k});
end
fprintf(fs, '%s\n', headers{end});

fprintf(fs, '%s\t', 'Sigma');
for index = 1:nbif-1,
    fprintf(fs, '%f\t', allResults{index}.p(1));
end
fprintf(fs, '%f\n', allResults{nbif}.p(1));

fprintf(fs, '%s\t', 'b');
for index = 1:nbif-1,
    fprintf(fs, '%f\t', allResults{index}.p(2));
end
fprintf(fs, '%f\n', allResults{nbif}.p(2));

fprintf(fs, '%s\t', 'a before');
for index = 1:nbif-1,
    fprintf(fs, '%f\t', allResults{index}.p(3));
end
fprintf(fs, '%f\n', allResults{nbif}.p(3));

fprintf(fs, '%s\t', 'a middle');
for index = 1:nbif-1,
    fprintf(fs, '%f\t', allResults{index}.p(4));
end
fprintf(fs, '%f\n', allResults{nbif}.p(4));

fprintf(fs, '%s\t', 'a after');
for index = 1:nbif-1,
    fprintf(fs, '%f\t', allResults{index}.p(5));
end
fprintf(fs, '%f\n', allResults{nbif}.p(5));

fclose(fs);

end
