%Preprocess RNAseq data using the following steps
%1. log2 tranform
%2. Filter out low-expressed genes
%3. Select highly variable geness
%4. if needed, infer cell stage by using prinicipal curve analysis.

function RNAseq_preprocess(dataset, log_mode, pseudotime_mode);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

lowgene_threshold = 1; %threshold value for low expressed genes
lowgene_fraction_max = 0.7; % maximum fraction of low-expressed cells allowed for each gene
ngene_select = 1000; %number of highly variable genes selected

if ~exist('log_mode'),
    log_mode = 1;
end

if ~exist('pseudotime_mode'),
    pseudotime_mode = 0;
end

txtout_mode = 1;

fs = fopen(dataFile);
header1 = fgetl(fs);
I = findstr(header1, char(9));
ncell = length(I);
cell_id = cell(1, ncell);
for k = 1:ncell-1,
    cell_id{k} = header1(I(k)+1:I(k+1)-1);
end
cell_id{end} = header1(I(end)+1:end);

if pseudotime_mode == 0, %if the temporal information is given, this should be included in the second row
    header2 = fgetl(fs);
    I = findstr(header2, char(9));
    cell_stage = zeros(1, ncell);
    for k = 1:ncell-1,
        cell_stage(k) = str2num(header2(I(k)+1:I(k+1)-1));
    end
    cell_stage(end) = str2num(header2(I(end)+1:end));
end

fmt = '%s';
for k = 1:ncell,
    fmt = [fmt, '%f'];
end
c = textscan(fs, fmt, 'delimiter', char(9));

gname = c{1};
clear pro;
pro.gname = gname;
ngene = length(gname);
pro.cell_stage = cell_stage;
pro.cell = cell_id;
pro.expr = zeros(ncell, ngene);
for k = 1:ncell,
    pro.expr(k, :) = c{k+1}';
end

islowgene = mean(pro.expr<lowgene_threshold,1) > lowgene_fraction_max;
pro.expr = pro.expr(:, ~islowgene);
pro.gname = pro.gname(~islowgene);
fanoFactors = var(pro.expr)./mean(pro.expr);
[sortedFanoFactors, sortedFanoFactorsIndices] = sort(fanoFactors,'descend');
pro.expr = pro.expr(:,sortedFanoFactorsIndices(1:ngene_select));
pro.gname = pro.gname(sortedFanoFactorsIndices(1:ngene_select));

if log_mode == 1,
    pro.expr = log2(pro.expr+1);
end

save(processDataMat, 'pro');

if pseudotime_mode == 1,
    EstimatePseudotime(dataset);
end

if txtout_mode == 1,
    %output data
    fout = fopen(processDataTxt, 'w+');
    fprintf(fout, '%s\n', header1);
    
    fprintf(fout, '%s\t', 'Stage');
    for k = 1:ncell-1
        fprintf(fout, '%f\t', pro.cell_stage(k));
    end
    fprintf(fout, '%f\n', pro.cell_stage(end));
    
    for k = 1:ngene_select,
        fprintf(fout, '%s\t', pro.gname{k});
        for j = 1:ncell-1,
            fprintf(fout, '%f\t', pro.expr(j, k));
        end
        fprintf(fout, '%f\n', pro.expr(ncell, k));
    end
    
    fclose all;
end

end
