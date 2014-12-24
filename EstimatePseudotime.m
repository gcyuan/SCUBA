function EstimatePseudotime(dataset, pcvmethod, anchorGene)
% options for pcvmethod are 'Rprincurve' and 'ksegments'
% use anchorGene to order pseudotime so that the average value of 
% anchorGene in the first 1000 cells is higher than for the last 1000 cells

% set default pcvmethod
if ~exist('pcvmethod', 'var')
    pcvmethod = 'Rprincurve';
end

% set default anchorGene
if ~exist('anchorGene', 'var')
    anchorGene = false;
end

[~, processDataMat, ~, ~, ~, ~, intermediate_filesDir, ~] = initialization(dataset);

load(processDataMat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1: dimension reduction by using tSNE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndim = 3;
method = 'tSNE';
[mappedX, ~] = compute_mapping(pro.expr, method, ndim);
pro.tsne = mappedX;
save(processDataMat, 'pro');

if strcmp(pcvmethod, 'Rprincurve')
    
    % save csv file to process with R's princurve
    ncell = size(pro.expr, 1);
    csvout = fullfile(intermediate_filesDir, [dataset '_tsne_d' num2str(ndim) '.csv']);
    fout = fopen(csvout, 'w+');
    ntsne = size(pro.tsne, 2);
    for k = 1:ntsne-1,
        fprintf(fout, '%s,', ['T', int2str(k)]);
    end
    fprintf(fout, '%s\n', ['T', int2str(ntsne)]);
    for k = 1:ncell,
        for j = 1:ntsne-1,
            fprintf(fout, '%f,', pro.tsne(k, j));
        end
        fprintf(fout, '%f\n', pro.tsne(k, end));
    end
    fclose(fout);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2: Fit tSNE data by Principal curve analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch pcvmethod
    case 'Rprincurve'
        % Calling R's principal curve
        system(['Rscript principal_curve_analysis_tsne_d3.R ' dataset]);
        
        % pcvin = fullfile(intermediate_filesDir, [dataset '_tsne_d' num2str(ndim) '_pcv.csv']);
        parin = fullfile(intermediate_filesDir, [dataset '_tsne_d' num2str(ndim) '_lambda.csv']);
        fin3 = fopen(parin);
        s = fgetl(fin3);
        c3 = textscan(fin3, '%s%f', 'delimiter', ',');
        fclose(fin3);
        pro.pseudotime = c3{2};
    case 'ksegments'       
        k_max = 10;
        alpha = 2;
        lambda = 1;
        INT_PLOT = 1;
        [edges,vertices,of,y]=k_seg_soft(pro.tsne,k_max,alpha,lambda,INT_PLOT);
        pro.pseudotime = y(:, 1);
end

if anchorGene
    % We use anchorGene to find beginning of pseudotime
    [~, loc] = ismember(anchorGene, pro.gname);
    [~, pseudotimeSortIndices] = sort(pro.pseudotime);
    nCellsAnchor = 1000;
    anchorGeneAverageBeginning = mean(pro.expr(pseudotimeSortIndices(1:nCellsAnchor),loc));
    anchorGeneAverageEnd = mean(pro.expr(pseudotimeSortIndices(end-nCellsAnchor+1:end),loc));

    if anchorGeneAverageEnd > anchorGeneAverageBeginning
        % We flip time
        pro.pseudotime = max(pro.pseudotime) - pro.pseudotime;
    end
end

tmin = min(pro.pseudotime);
tmax = max(pro.pseudotime);
tbin = 8;
pro.cell_stage = ceil(tbin*(pro.pseudotime - tmin + 0.0001)/(tmax + 0.0002 - tmin));

save(processDataMat, 'pro');