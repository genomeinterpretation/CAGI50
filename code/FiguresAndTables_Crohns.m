%% Generate the figures and tables for Crohns after running generateResults_Crohns.m
%% The table is stored in '../results' folder.

addpathscript;
s = load(['../results/Crohns.mat']);
Results = s.Results;
DSAbbrs = fieldnames(s.Results);
tableName = 'Crohns';
dsAbbr = DSAbbrs{1};
ds = Results.(dsAbbr).ds;
perfs = Results.(dsAbbr).perfs;
n_main_methods = 2;

selectedPerfs = selectMethodsAndBaselines(perfs, ds, n_main_methods);


figure;
rocType = 'standard';
ROCPlot(selectedPerfs, ds, rocType, []);
disp('Standard ROC plot created.')

figure;
rocType = 'standard';
densityPlot(selectedPerfs(1), ds);
disp('Density plot created.')

figure;
show_llrp = true;
smoothing = true;
priors = [];
posteriorAndllrPlusPlot(selectedPerfs(1), ds, priors, show_llrp, smoothing)
disp('local lr+ plot created.')

figure;
smoothing = true;
prior = 0.013;
RRPlot(selectedPerfs(1), ds, prior, smoothing)
disp('Posterior and local lr+ plot created.')

figure;


priors = selectedPerfs(1).priors;


generateTable(Results, tableName);
disp('Table created in ../results');
