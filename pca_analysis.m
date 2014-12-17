%Process gene expression data by PCA analysis.

function pca_analysis(dataset);

[dataFile processDataMat processDataTxt PCAdataFile dataFolder resultsDir intermediate_filesDir figuresDir] = initialization(dataset);

load(processDataMat);

%PCA analysis using all samples
ncell = size(pro.expr, 1);
ngene = size(pro.expr, 2);
npca = min(ncell-1,ngene); %meaningful number of PCs

mu = mean(pro.expr, 1);
expr = pro.expr - repmat(mu, ncell, 1);
[c, s] = pca(expr);
expr_all = pro.expr - repmat(mu, ncell, 1);
s_all = expr_all * c;
X = s_all(:, 1:npca); %reduced data
pro.pca = X;
pro.weight = c;

%PCA analysis using samples from the last cell stage.
I = find(pro.cell_stage == max(pro.cell_stage)); 
npca2 = min(length(I)-1,ngene); %meaningful number of PCs

mu = mean(pro.expr(I, :), 1);
nI = length(I);
expr = pro.expr(I,:) - repmat(mu, nI, 1);
[c, s] = pca(expr); 
expr_all = pro.expr - repmat(mu, ncell, 1);
s_all = expr_all * c;
X = s_all(:, 1:npca2); %reduced data
pro.pca2 = X;
pro.weight2 = c;

save(PCAdataFile, 'pro');

end