function MassCytometry_preprocess(dataset, select_marker_names, pcvmethod, anchorGene, seed)
%Preprocess Mass-cytometry data using the following steps
%1. convert data from fcs to mat format
%2. infer cell stage by using prinicipal curve analysis.
% options for pcvmethod are 'Rprincurve' (default) and 'ksegments'

% option use a seed to ensure reproducibility of the results
if exist('seed', 'var')
    rng(seed)
end

% set default pcvmethod
if ~exist('pcvmethod', 'var')
    pcvmethod = 'Rprincurve';
end

% set default anchorGene
if ~exist('anchorGene', 'var')
    anchorGene = false;
end

[~, processDataMat, ~, ~, dataFolder, ~, ~, ~] = initialization(dataset);

fcsfile = fullfile(dataFolder, [dataset, 'Data.fcs']);

[fcsdat, fcshdr, ~] = fca_readfcs(fcsfile);
npar = length(fcshdr.par);
select_col = 1:npar;
nmarker = length(select_col);
marker_names = cell(1, nmarker);
for k = 1:nmarker,
    marker_names{k} = fcshdr.par(select_col(k)).name;
end

M = fcsdat(:, select_col);

ncell = size(M,1);

cell_id = cell(ncell,1);
for index = 1:ncell
    cell_id{index,1} = ['cell_' num2str(index)];
end

pro.cell = cell_id;

if nargin == 1
    % We select all markers
    select_marker = 3:(size(M,2)-1); %remove non-gene related parameters, such as time and cell_length.
else
    % We select a subset of markers
    [~, select_marker] = ismember(select_marker_names, marker_names);
end

pro.expr = M(:, select_marker);
pro.gname = marker_names(select_marker);

if exist('seed', 'var')
    pro.seed = seed;
end
pro.pcvmethod = pcvmethod;
pro.anchorGene = anchorGene;

save(processDataMat, 'pro');

%Estimate the pseudotime
EstimatePseudotime(dataset, pcvmethod, anchorGene);