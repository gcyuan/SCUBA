function getBifurcationDirections
% Date 5/14/14: written by Eugenio Marco
% Calculate and compare bifurcation directions for deng2014 and guo2010 datasets

[~, ~, resultsDir, intermediate_filesDir, ~] = initialization('deng2014');

tree_in = fullfile(intermediate_filesDir,'mixture_tree.mat');
load(tree_in)

clusters_in = fullfile(intermediate_filesDir,'mixture_cluster.mat');
load(clusters_in)

genefile = fullfile(intermediate_filesDir,'gene_expr.mat');
load(genefile);


bifurcationParents = T.pa(T.pa(1:end-1)==T.pa(2:end));

% bifurcationDirections = zeros(length(bifurcationParents), size(pro.expr,2));
bifurcationDirectionsUni = zeros(length(bifurcationParents), size(pro.expr,2));
for index = 1:length(bifurcationParents)
    bifurcationChildren = T.clu_id(T.pa == bifurcationParents(index));
    cellsCluster1 = y.clu_idx == bifurcationChildren(1);
    mu1 = mean(pro.expr(cellsCluster1,:));
    cellsCluster2 = y.clu_idx == bifurcationChildren(2);
    mu2 = mean(pro.expr(cellsCluster2,:));
    vuni = mu1 - mu2;
    vuni = vuni/norm(vuni);
    bifurcationDirectionsUni(index,:) = vuni;
end

% We flip bifurcation 2 for deng2014 to agree with the order in guo2010
% bifurcationDirectionsUni(2,:) = -bifurcationDirectionsUni(2,:);


genes = pro.gname;

%%
figure
plot(bifurcationDirectionsUni(1,:),bifurcationDirectionsUni(2,:),'.')
set(gca,'fontsize',18)
xlabel('X32')
ylabel('X64')
xlim([-.17 .17])
ylim([-.2 .1])
title('Gene Weights for RNA-seq Data')
%%
% Read Guoji's data
perturbationToPCA48 = myReadPerturbationToPCA;

bifurcationDirections48 = zeros(2,48);
allIndices = [];
% Subset data
for index = 1:length(perturbationToPCA48.gene)
    [~, ib, ~] = intersect(genes,perturbationToPCA48.gene(index));
    if ~isempty(ib)
        bifurcationDirections48(:,index) = bifurcationDirectionsUni(:,ib);
        allIndices = [allIndices, index];
    else
        bifurcationDirections48(:,index) = 0;
    end
end
bifurcationDirections48 = bifurcationDirections48(:,allIndices);

fs=10;
%% Like in Fig 2b
figure
plot(perturbationToPCA48.X32,perturbationToPCA48.X64,'o')
set(gca,'fontsize',18)

for index = 1:length(perturbationToPCA48.X32)
    text((.0252+perturbationToPCA48.X32(index)),...
        (.0252+perturbationToPCA48.X64(index)),...
        perturbationToPCA48.gene{index},'fontsize',fs,'color','k')
end
title('Gene Weights for RT-PCR Data')

xlabel('X32')
ylabel('X64')
xlim([-.4 .4])
ylim([-.4 .4])
%%
y32 = perturbationToPCA48.X32(allIndices);
x32 = bifurcationDirections48(1,:)';
[r32, p32] = corrcoef(x32,y32);
disp(['32C: Pearson''s correlation = ' num2str(r32(1,2))])

[b32,bint32,res32,rint32,stats32] = regress(y32, [ones(size(allIndices))'  x32]);
%%
figure
plot(bifurcationDirections48(1,:),perturbationToPCA48.X32(allIndices),'o')
xl = get(gca,'xlim');
yl = get(gca,'ylim');
hold on

for index = 1:length(x32)
    text((.0052+x32(index)),...
        (.0052+y32(index)),...
        perturbationToPCA48.gene{allIndices(index)},'fontsize',16,'color','k')    
end
xFine = xl(1):.01:xl(2);
plot(xFine,b32(1)+b32(2)*(xFine))
xlim(xl)
ylim(yl)
set(gca,'fontsize',18)
xlabel('32-Cell Bifurcation from RNA-Seq')
ylabel('32-Cell Bifurcation from RT-PCR')
title(['\it{R}^2 = ' num2str(stats32(1),2)])


filenameFig = fullfile(pwd,'figures','figCorr32.eps');

hgexport(gcf, filenameFig)
%%
y64 = perturbationToPCA48.X64(allIndices);
x64 = bifurcationDirections48(2,:)';
[r64, p64] = corrcoef(x64,y64);
disp(['64C: Pearson''s correlation = ' num2str(r64(1,2))])

% cov(bifurcationDirections48(1,allIndices),perturbationToPCA48.X32(allIndices)')
[b64,bint64,res64,rint64,stats64] = regress(y64, [ones(size(allIndices))'...
    x64]);
figure
plot(x64, y64 ,'o')
xl = get(gca,'xlim');
yl = get(gca,'ylim');
hold on

for index = 1:length(x64)
    text((.0052+x64(index)),...
        (.0052+y64(index)),...
        perturbationToPCA48.gene{allIndices(index)},'fontsize',16,'color','k')
end
xx = xl(1):.01:xl(2);
plot(xx,b64(1)+b64(2)*xx)
% xlim([-0.06 0.04])
xlim(xl)
ylim(yl)
set(gca,'fontsize',18)
xlabel('64-Cell Bifurcation from RNA-Seq')
ylabel('64-Cell Bifurcation from RT-PCR')
title(['\it{R}^2 = ' num2str(stats64(1),2)])

% title(['64-Cell Correlation = ' num2str(r64(1,2),3)])

% cov(bifurcationDirections48(2,allIndices),perturbationToPCA48.X64(allIndices)')

filenameFig = fullfile(pwd,'figures','figCorr64.eps');

hgexport(gcf, filenameFig)

%% Print ranked bifurcations for RNA-seq data
file1 = fullfile(resultsDir,'ranked_bifurcation1.rnk');
file2 = fullfile(resultsDir,'ranked_bifurcation2.rnk');

fs1 = fopen(file1, 'W+');
fprintf(fs1, '# Ranked list for first bifurcation\n');
for index = 1:length(genes)
    fprintf(fs1, '%s\t', genes{index});
    fprintf(fs1, '%d\n', bifurcationDirectionsUni(1,index));
end
fclose(fs1);

fs2 = fopen(file2, 'W+');
fprintf(fs2, '# Ranked list for second bifurcation\n');
for index = 1:length(genes)
    fprintf(fs2, '%s\t', genes{index});
    fprintf(fs2, '%d\n', bifurcationDirectionsUni(2,index));
end
fclose(fs2);

%%

file3 = fullfile(resultsDir,'RNA-seq_gene_weights_bifurcations.txt');

fs1 = fopen(file3, 'W+');
fprintf(fs1, 'Gene\tX32\tX64\n');
for index = 1:length(genes)
    fprintf(fs1, '%s\t', genes{index});
    fprintf(fs1, '%d\t', bifurcationDirectionsUni(1,index));
    fprintf(fs1, '%d\n', bifurcationDirectionsUni(2,index));
end
fclose(fs1);

%%
file4 = fullfile(resultsDir,'RT-PCR_gene_weights_bifurcations.txt');

fs1 = fopen(file4, 'W+');
fprintf(fs1, 'Gene\tX32\tX64\n');
for index = 1:length(perturbationToPCA48.X32)
    fprintf(fs1, '%s\t', perturbationToPCA48.gene{index});
    fprintf(fs1, '%d\t', perturbationToPCA48.X32(index));
    fprintf(fs1, '%d\n', perturbationToPCA48.X64(index));
end
fclose(fs1);
%% Print ranked bifurcations for RT-PCR data
[~, ~, resultsDirGuo, ~, ~] = initialization('guo2010');

file1 = fullfile(resultsDirGuo,'ranked_bifurcation1.rnk');
file2 = fullfile(resultsDirGuo,'ranked_bifurcation2.rnk');

fs1 = fopen(file1, 'W+');
fprintf(fs1, '# Ranked list for first bifurcation\n');
for index = 1:length(perturbationToPCA48.X32)
    fprintf(fs1, '%s\t', perturbationToPCA48.gene{index});
    fprintf(fs1, '%d\n', perturbationToPCA48.X32(index));
end
fclose(fs1);

fs2 = fopen(file2, 'W+');
fprintf(fs2, '# Ranked list for second bifurcation\n');
for index = 1:length(perturbationToPCA48.X32)
    fprintf(fs2, '%s\t', perturbationToPCA48.gene{index});
    fprintf(fs2, '%d\n', perturbationToPCA48.X64(index));
end
fclose(fs2);
%%
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFile, dataFolder, resultsDir, intermediate_filesDir, figuresDir] = initialization(dataset)

% set path of data files
dataFile = [dataset 'Data.txt'];
dataFolder = fullfile(pwd,'sample_data',dataset);
save('dataFolder.mat','dataFolder');

% Alternatively, uncomment the following line and select interactively
% [dataFile dataFolder] = setDataFile;

% Create results folder in data folder if not previously present
% createSubfolders(dataFolder)

% Set the results folder variable
resultsDir = fullfile(dataFolder,'results');
intermediate_filesDir = fullfile(dataFolder,'intermediate_files');
figuresDir = fullfile(dataFolder,'figures');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function perturbationToPCA = myReadPerturbationToPCA
% Date 5/14/14: written by Eugenio Marco
% Read gene perturbation effects on bifurcation axis

[dataFile, dataFolder, resultsDir, intermediate_filesDir, figuresDir] = initialization('guo2010');

txtfile = fullfile(intermediate_filesDir, 'perturbation_to_mixture.txt');

fs = fopen(txtfile);
s = fgetl(fs);
I = findstr(s, char(9));
nColumns = length(I);

fmt = '%s';
for k = 1:nColumns,
    fmt = [fmt, '%f'];
end

data = textscan(fs, fmt, 'delimiter', char(9));

perturbationToPCA.gene = data{1};
perturbationToPCA.X32 = data{2};
perturbationToPCA.X64 = data{3};
perturbationToPCA.PCAs = cell2mat(data(4:end));
fclose(fs);

end