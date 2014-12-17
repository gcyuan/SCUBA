function EstimatePseudotime(dataset);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

load(processDataMat);

%step 1: dimension reduction by using tSNE
ndim = 3;
method = 'tSNE';
[mappedX, mapping] = compute_mapping(pro.expr, method, ndim);
pro.tsne = mappedX;

save(processDataMat, 'pro');

%step 2: Fit tSNE data by Principal curve analysis
k_max = 10;
alpha = 2;
lambda = 1;
INT_PLOT = 1;
[edges,vertices,of,y]=k_seg_soft(pro.tsne,k_max,alpha,lambda,INT_PLOT);
pro.pseudotime = y(:, 1);
tmin = min(pro.pseudotime);
tmax = max(pro.pseudotime);
tbin = 8;
pro.cell_stage = ceil(tbin*(pro.pseudotime - tmin + 0.0001)/(tmax + 0.0002 - tmin));

save(processDataMat, 'pro');