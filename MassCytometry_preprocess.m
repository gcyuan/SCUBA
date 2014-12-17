%Preprocess Mass-cytometry data using the following steps
%1. convert data from fcs to mat format
%2. infer cell stage by using prinicipal curve analysis.

function MassCytometry_preprocess(dataset);

rng(7)

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

fcsfile = fullfile(dataFolder, [dataset, 'Data.fcs']);

[fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(fcsfile);
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

select_marker = 3:(size(M,2)-1); %remove non-gene related parameters, such as time and cell_length.
pro.expr = M(:, select_marker);
pro.gname = marker_names(select_marker);

save(processDataMat, 'pro');

%Estimate the pseudotime
EstimatePseudotime(dataset);