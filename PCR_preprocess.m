%Preprocessing for PCR data.
%1. data converted to matlab format
%2. Filter out genes whose expression is undetectable in most cells.
%3. If needed, infer cell stage by using principal curve analysis.

function PCR_preprocess(dataset, pseudotime_mode);

if ~exist('pseudotime_mode'),
    pseudotime_mode = 0;
end

txtout_mode = 1; %set to 1 if also want to output the data in tab-delimited format. 

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

lowgene_fraction_max = 0.8; % maximum fraction of low-expressed cells allowed for each gene

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
if exist(cell_stage),
    pro.cell_stage = cell_stage;
end
pro.cell = cell_id;
pro.expr = zeros(ncell, ngene);
for k = 1:ncell,
    pro.expr(k, :) = c{k+1}';
end

islowgene = mean(pro.expr==0,1) > lowgene_fraction_max;
pro.expr = pro.expr(:, ~islowgene);
pro.gname = pro.gname(~islowgene);

save(processDataMat, 'pro');

if pseudotime_mode == 1,
    EstimatePseudotime(dataset);
end

if txtout_mode == 1,
    if pseudotime_mode == 1,
        load(processDataMat);
    end
    fout = fopen(processDataTxt, 'w+');
    fprintf(fout, '%s\n', header1);
    
    fprintf(fout, '%s\t', 'Stage');
    for k = 1:ncell-1
        fprintf(fout, '%f\t', pro.cell_stage(k));
    end
    fprintf(fout, '%f\n', pro.cell_stage(end));
    
    ngene_select = length(pro.gname);
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
