%Initial guess of cell hierarchy, using k-means clustering.

function initialize_tree(dataset, cluster_mode);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

load(PCAdataFile);

outfile = fullfile(intermediate_filesDir,'initial_tree.mat');

clear y;
y.cell = pro.cell;
y.gname = pro.gname;
y.cell_stage = pro.cell_stage;
y.expr = pro.expr;  %expression level
y.pca = pro.pca;    %pca coordinate
y.pca2 = pro.pca2;  %pca coordinate based on the last time step

ncell = size(y.expr, 1);
stage = sort(unique(y.cell_stage));
nstage = length(stage);

y.cm = cell(nstage, 1);  %cluster mean
y.cpt = cell(nstage, 1); %parental cluster
y.clu_id = cell(nstage, 1); %Clusters that are associated with each stage
y.nclu = zeros(nstage, 1);  %number of clusters at each stage

clear pro;

%Preparation for computing gap-statistics
if 1,
    nsamp = 1000;
    switch cluster_mode,
        case 'pca',
            ndim = size(y.pca, 2);
            r = rand(nsamp, ndim);
            w = diag(std(y.pca));
            r = r*w/sqrt(1/12);
        case 'pca2',
            ndim = size(y.pca2, 2);
            r = rand(nsamp, ndim);
            w = diag(std(y.pca2));
            r = r*w/sqrt(1/12);
        otherwise,
            ndim = size(y.expr, 2);
            r = rand(nsamp, ndim);
            w = diag(std(y.expr));
            r = r*w/sqrt(1/12);
    end
    
    kmax = 10;
    Wk0 = zeros(1, kmax);
    for k = 1:kmax,
        [idx, c, sumd] = kmeans(r, k, 'replicates', 20);
        Wk0(k) = sum(sumd)/nsamp;
    end
end

t0 = 1;
I = find(y.cell_stage == stage(t0));
switch cluster_mode,
    case 'pca',
        X = y.pca(I, :);
    case 'pca2',
        X = y.pca2(I, :);
    otherwise,
        X = y.expr(I, :);
end
N = length(I);

%Estimate the number of clusters at initial stage.
kmax = min(10, length(I));
Wk = zeros(1, kmax);
G = zeros(1, kmax);
for k = 1:kmax,
    [idx, c, sumd] = kmeans(X, k, 'replicates', 20);
    Wk(k) = sum(sumd)/nsamp;
    G(k) = Wk0(k) - Wk(k);
end
[Gmax, imax] = max(G);
[idx, c, sumd] = kmeans(X, imax, 'replicates', 20);
    
y.cm{t0} = c;
y.cpt{t0} = 0;
y.clu_id{t0} = 1:imax;
y.nclu(t0) = imax;
y.clu_idx(I) = idx;
clu_cnt = imax;

%Define a few parameters in the code
min_cell_number = 15; %Minimum size of a cluster to be split

min_percentage_cells_split = 0.25; %Minimum fraction of cells in the smaller cluster during a bifurcation

Wk = zeros(nstage, 5);

for t = t0+1:nstage,
    t
    I = find(y.cell_stage == stage(t));
    switch cluster_mode,
        case 'pca',
            X = y.pca(I, :);
        case 'pca2',
            X = y.pca2(I, :);
        otherwise,
            X = y.expr(I, :);
    end
    
    N = size(X, 1);
    nclu_up = y.nclu(t-1); %Number of clusters at previous time step
    
    %map cells to correspond clusters
    if nclu_up == 1,
        clumap = ones(N, 1);
    else,
        clumap = zeros(N, 1);
        D = zeros(N, nclu_up);
        for k = 1:nclu_up,
            r = X - repmat(y.cm{t-1}(k, :), N, 1);
            D(:, k) = sqrt(sum(r.*r, 2));
        end
        [Dmin, Imin] = min(D,[], 2);
        clumap = Imin;
    end
    
    % test bifurcation
    y.nclu(t) = 0;
    y.cm{t} = [];
    y.clu_id{t} = [];
    y.cpt{t} = [];
    
    for j = 1:nclu_up,
        L = find(clumap == j);
        if isempty(L),
            warning(['Empty Cluster mapping encountered at t = ', int2str(t)]);
            continue;
        else
            XX = X(L, :);
            if length(L) < min_cell_number, %split only if there are enough cells
                y.nclu(t) = y.nclu(t) + 1;
                y.clu_id{t} = [y.clu_id{t} clu_cnt+1];
                y.cm{t} = [y.cm{t}; mean(XX, 1)];
                y.cpt{t} = [y.cpt{t} y.clu_id{t-1}(j)];
                y.clu_idx(I(L)) = clu_cnt+1;
                clu_cnt = clu_cnt + 1;
            else
                Wk = zeros(1, 2);
                G = zeros(1, 2);
                idx = cell(1, 2);
                c = cell(1, 2);
                for k = 1:2,
                    [idx{k}, c{k}, sumd] = kmeans(XX, k, 'replicates', 20); %sumd is the sum of square.
                    Wk(k) = sum(sumd)/length(L);
                    G(k) = Wk0(k) - Wk(k);
                end
                if G(2) > G(1)&& min([sum(idx{2}==1) sum(idx{2}==2)])/length(idx{2}) >  min_percentage_cells_split
                    y.nclu(t) = y.nclu(t) + 2;
                    y.clu_id{t} = [y.clu_id{t} clu_cnt+(1:2)];
                    y.cm{t} = [y.cm{t}; c{2}];
                    y.cpt{t} = [y.cpt{t} y.clu_id{t-1}(j)*ones(1,2)];
                    y.clu_idx(I(L)) = clu_cnt + idx{2};
                    clu_cnt = clu_cnt + 2;
                else
                    y.nclu(t) = y.nclu(t) + 1;
                    y.clu_id{t} = [y.clu_id{t} clu_cnt+1];
                    y.cm{t} = [y.cm{t}; c{1}];
                    y.cpt{t} = [y.cpt{t} y.clu_id{t-1}(j)];
                    y.clu_idx(I(L)) = clu_cnt + idx{1};
                    clu_cnt = clu_cnt + 1;
                end
            end
        end
    end
end

save(outfile, 'y');



