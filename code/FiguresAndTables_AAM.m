%% Generate the figures and tables for all Annotate all missense datasets.
%% after running generateResults_AAM.m
%% The table is stored in '../results' folder.
% Because HGMD variants are not provided, the script only generates the figures 
% and tables corresponding to the AAM1CV and AAM2CV sheets in Table 4. The 
% figures and tables for the bi-class genes generated from this code does not
% correspond to any bi-class sheet in Table 4 because of the missing HGMD variants.
addpathscript;
s = load(['../results/AAM.mat']);
DS = fieldnames(s.Results);
tableName = 'AAM';
Results = struct();
for d = 1:length(DS)
    dsAbbr = DS{d};
    Results.(dsAbbr) = s.Results.(dsAbbr);
    perfs = Results.(dsAbbr).perfs;
    ds = Results.(dsAbbr).ds;
    n_main_methods = 2;
    selectedPerfs = selectMethodsAndBaselines(perfs, ds, n_main_methods); 
    

    figure;
    rocType = 'standard';
    ROCPlot(selectedPerfs, ds, rocType, []);
    disp('Standard ROC plot created.')

    figure;
    rocType = 'truncated';
    priors = 0.1;
    ROCPlot(selectedPerfs, ds, rocType, priors); 
    disp('Truncated ROC plot created.')

    figure;
    rocType = 'log-log';
    ROCPlot(selectedPerfs, ds, rocType, priors); 
    disp('log-log ROC plot created.')

    figure;
    show_llrp = true;
    smoothing = true;
    priors = selectedPerfs(1).priors;
    posteriorAndllrPlusPlot(selectedPerfs(1), ds, priors, show_llrp, smoothing)
    disp('Posterior and local lr+ plot created.')

end
generateTable(Results, tableName);
disp('Table created in ../results');
