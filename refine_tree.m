%Refine tree structure by maximizing penalized likelihood method

function refine_tree(dataset, cluster_mode);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

infile = fullfile(intermediate_filesDir,'initial_tree.mat');
outfile = fullfile(intermediate_filesDir,'final_tree.mat');
txtoutfile = fullfile(intermediate_filesDir,'final_tree.txt');

load(infile);
stage = sort(unique(y.cell_stage));
nstage = length(stage);

clear T %this is the main tree

%initialze the tree
nclu_total = sum(y.nclu);

%initialize time, cluster_id, and parent_id
T.t = zeros(nclu_total, 1);
T.stage = zeros(nclu_total, 1);
T.clu_id = zeros(nclu_total, 1);
T.pa = zeros(nclu_total, 1); %parental clusters
tstart = 1;
for k = 1:nstage,
    nclu = y.nclu(k);
    tend = tstart + nclu - 1;
    T.t(tstart:tend) = k;
    T.stage(tstart:tend) = stage(k);
    T.clu_id(tstart:tend) = y.clu_id{k};
    T.pa(tstart:tend) = y.cpt{k};
    tstart = tend + 1;
end

switch cluster_mode,
    case 'pca',
        X = y.pca;
    case 'pca2',
        X = y.pca2;
    otherwise,
        X = y.expr;
end

%initialize mu: position of the node.
ndim = size(X, 2); %dimension of the problem
T.mu = zeros(nclu_total, ndim);
tstart = 1;
for k = 1:nstage,
    nclu = y.nclu(k);
    tend = tstart + nclu - 1;
    T.mu(tstart:tend, :) = y.cm{k};
    tstart = tend + 1;
end

%initilalize sigma: variance within each cluster
%%%Should we use the variance of each cluster as initialization?
S = ones(1, ndim);
T.sigma = repmat(S, nclu_total, 1); %assume that the sigma is constant across clusters

%initilalize D: temporal continuity of the clusters
lambda = 0.5; %diffusion constant for cluster drift.


%initialize cluster assignment for each cell.
T.s = y.clu_idx;


%Update tree structure based on penalized likelihood
%Main iteration

iter_stop = 0;
tol = 0.000000001; %criterion for convergence
maxiter = 10; %maximum number of iterations
iter = 1;
while iter < maxiter & iter_stop == 0,
    %update s
    s_old = T.s;
    for k = 1:nstage,
        clu_idx = find(T.t == k);
        mu = T.mu(clu_idx, :);
        sigma = T.sigma(clu_idx, :);
        nclu = length(clu_idx);
        if nclu > 1, %only update if there are more than one clusters
            I = find(y.cell_stage == stage(k));
            nI = length(I);
            d = zeros(nI, nclu);
            for j = 1:nclu,
                %                     M = (y.pca(I, :) - repmat(mu(j, :), nI, 1))./repmat(sigma(j, :), nI, 1);
                %                     M = (y.expr(I, :) - repmat(mu(j, :), nI, 1))./repmat(sigma(j, :), nI, 1);
                M = (X(I, :) - repmat(mu(j, :), nI, 1))./repmat(sigma(j, :), nI, 1);
                d(:, j) = sqrt(sum(M.*M, 2));
            end
            [dmin, idx] = min(d, [], 2);
            T.s(I) = clu_idx(idx);
            
            if unique(idx) < nclu, %make adjustment if empty cluster emerges
                warning(['Empty clusters identified at t = ', int2str(k)]);
                nidx = zeros(1, nclu);
                for j = 1:nclu,
                    nidx(j) = length(find(idx == j));
                end
                J = find(nidx == 0);
                for j = 1:length(J), %make sure no empty cluster exists by adding the closest point to the cluster
                    [zmin, zidx] = min(d(:, J(j)));
                    T.s(zidx) = clu_idx(J(j));
                end
            end
        end
    end
    
    %update mu
    mu_old = T.mu;
    for j = 1:nclu_total,
        t = T.t(j);
        J = find(T.s == j);
        nJ = length(J);
        % %             sumx = sum(y.pca(J, :), 1);
        %             sumx = sum(y.expr(J, :), 1);
        sumx = sum(X(J, :), 1);
        j0 = T.pa(j);
        if j0 == 0,
            T.mu(j, :) = sumx./nJ;
        else,
            mu0 = T.mu(j0, :);
            T.mu(j, :) = (2*lambda*mu0 + sumx)./(2*lambda + nJ);
        end
    end
    
    
    %update pa:
    
    pa_old = T.pa;
    for k = 2:nstage,
        clu_idx = find(T.t == k);
        nclu = length(clu_idx);
        mu = T.mu(clu_idx, :);
        sigma = T.sigma(clu_idx, :);
        
        clu_idx_up = find(T.t == k-1);
        nclu_up = length(clu_idx_up);
        mu_up = T.mu(clu_idx_up, :);
        
        
        d = zeros(nclu_up, nclu);
        for i = 1:nclu_up,
            for j = 1:nclu,
                d(i,j) = sum((mu_up(i, :) - mu(j, :)).^2);
            end
        end
        [dmin, Imin] = min(d, [], 1);
        T.pa(clu_idx) = clu_idx_up(Imin);
        
        %Test if the configuration is legal
        npa = zeros(1, nclu_up);
        for i = 1:nclu_up,
            npa(i) = length(find(T.pa == clu_idx_up(i)));
        end
        if max(npa > 2) | min(npa < 1),
            %%% Should we allow this constraint to be violated?
            warning(['Binary tree assumption violated at t = ', int2str(t), '. Use old configuration instead.']);
            T.pa(clu_idx) = pa_old(clu_idx);
        end
    end
    
    disp(['iteration # ', int2str(iter)]);
    
    %check iteration success
    if 1,
        mu_tol = 0.000000000000001;
        delta_mu = sum(sum(abs(T.mu - mu_old)))
        delta_s  = sum(abs(T.s - s_old))
        delta_pa = sum(abs(T.pa - pa_old))
        
        if delta_mu < mu_tol & delta_s == 0 & delta_pa == 0,
            disp('Iteration successfully ended!');
            iter_stop = 1;
        end
    end
    
    iter = iter + 1;
end

if license('test', 'bioinformatics_toolbox')
    % Using biograph to plot the Tree
    plotTree(T, dataset)
else
    disp('Bioinformatics toolbox not present. Tree will not be plotted.')
end

save(outfile, 'T');

%Save the tree structure in txt format
fs = fopen(txtoutfile, 'w+');

headers = cell(1, ndim+4);
nheader = length(headers);
headers{1} = 'Cluster ID';
headers{2} = 'Time';
headers{3} = 'Cell Stage';
headers{4} = 'Parent cluster';

switch cluster_mode,
    case 'pca',
        for k = 1:ndim,
            headers{4+k} = ['MU PC ', int2str(k)];
        end
    case 'pca2',
        for k = 1:ndim,
            headers{4+k} = ['MU PC ', int2str(k)];
        end
    otherwise,
        for k = 1:ndim,
            headers{4+k} = ['MU ', y.gname{k}];
        end
end

for k = 1:nheader-1,
    fprintf(fs, '%s\t', headers{k});
end
fprintf(fs, '%s\n', headers{end});

for n = 1:nclu_total,
    fprintf(fs, '%d\t', n);
    fprintf(fs, '%d\t', T.t(n));
    fprintf(fs, '%d\t', T.stage(n));
    fprintf(fs, '%d\t', T.pa(n));
    for k = 1:ndim-1,
        fprintf(fs, '%e\t', T.mu(n, k));
    end
    fprintf(fs, '%e\n', T.mu(n, end));
end
fclose all;

end






