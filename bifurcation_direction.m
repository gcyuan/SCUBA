%Identify the bifurcation directions.
%Project data along the bifurcation directions.

function bifurcation_direction(dataset, cluster_mode);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

load(PCAdataFile);

tree_in = fullfile(intermediate_filesDir,'final_tree.mat');
bif_out = fullfile(intermediate_filesDir,'bifurcation_direction.mat');

DirectionOutfile = fullfile(intermediate_filesDir,'bifurcation_direction.txt');
ProjectionOutfile = fullfile(intermediate_filesDir,'projection_all_data.txt');

load(tree_in);

bifurcationEvents = find(T.pa(1:end-1) == T.pa(2:end));

if isempty(bifurcationEvents),
    warning('No bifurcation found.');
else,
    clear bif;
    bif.t = T.t(bifurcationEvents);
    bif.stage = T.stage(bifurcationEvents);
    ngene = length(pro.gname);
    ncell = length(pro.cell);
    nbif = length(bifurcationEvents);
         
    %projection of each cell onto the bifurcation direction
    bif.direction = zeros(ngene, nbif); %Bifurcation direction
    bif.proj = zeros(ncell, nbif);      %Projection of each cell along the bifurcation direction
    bif.clu_id = cell(1, nbif);         %corresponding cluster idx association with each bifurcation: [before, middle, after(1), after(2)] 
    for k = 1:nbif,
        I1 = bifurcationEvents(k);
        I2 = I1+1;
        pa = T.pa(I1);
        ppa = T.pa(pa);
        bif.clu_id{k} = [ppa, pa, I1, I2];
        t = T.t(I1);
        v = T.mu(I2, :) - T.mu(I1, :);
        switch cluster_mode,
            case 'pca',
                v = v/norm(v);
                bif.direction(:, k) = pro.weight*v';
                bif.proj(:, k) = pro.pca*v';
            case 'pca2'
                v = v/norm(v);
                bif.direction(:, k) = pro.weight2*v';
                bif.proj(:, k) = pro.pca2*v';
            otherwise,
                v = v/norm(v);
                bif.direction(:, k) = v;
                bif.proj(:, k) = pro.expr*v';
        end
    end
    
    save(bif_out, 'bif');
    
    %output the bifurcation direction.
    fs = fopen(DirectionOutfile, 'w+');
    
    fprintf(fs, '%s\t', 'Time');
    for k = 1:nbif-1,
        fprintf(fs, '%d\t', bif.t(k));
    end
    fprintf(fs, '%d\n', bif.t(end));
    
    fprintf(fs, '%s\t', 'Stage');
    for k = 1:nbif-1,
        fprintf(fs, '%d\t', bif.stage(k));
    end
    fprintf(fs, '%d\n', bif.stage(end));
    
    for j = 1:ngene,
        fprintf(fs, '%s\t', pro.gname{j});
        for k = 1:nbif-1,
            fprintf(fs, '%e\t', bif.direction(j, k));
        end
        fprintf(fs, '%e\n', bif.direction(j, end));
    end
    fclose(fs);
    
    %output the projection of data to bifurcation directions.
    fs = fopen(ProjectionOutfile, 'w+');
    
    fprintf(fs, '%s\t', 'Time');
    for k = 1:nbif-1,
        fprintf(fs, '%d\t', bif.t(k));
    end
    fprintf(fs, '%d\n', bif.t(end));
    
    fprintf(fs, '%s\t', 'Stage');
    for k = 1:nbif-1,
        fprintf(fs, '%d\t', bif.stage(k));
    end
    fprintf(fs, '%d\n', bif.stage(end));
    
    for j = 1:ncell,
        fprintf(fs, '%s\t', pro.cell{j});
        for k = 1:nbif-1,
            fprintf(fs, '%e\t', bif.proj(j, k));
        end
        fprintf(fs, '%e\n', bif.proj(j, end));
    end
    fclose all
end