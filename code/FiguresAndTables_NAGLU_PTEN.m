%% Generate the figures and tables for NAGLU and PTEN after running 
%% generateResults_NAGLU_PTEN.m.
%% The table is stored in '../results' folder.
addpathscript;
DS = {'NAGLU', 'PTEN'};
tableName = 'NAGLU_PTEN';
s = load(['../results/NAGLU_PTEN.mat']);
Results = struct();
for d = 1:length(DS)
    dsAbbr = DS{d};
    
    Results.(dsAbbr) = s.Results.(dsAbbr);
    perfs = Results.(dsAbbr).perfs;
    ds = Results.(dsAbbr).ds;
   
    topPerfs = selectMethodsAndBaselines(perfs, ds, 2); 
    
    
    figure;
    prior = perfs(1).data_prior;
    scatterPlot(topPerfs(1), ds, prior);
    disp('Scatter plot created.')

    figure;
    showKendall = true;
    showPearson = true;
    showSpearman = false;
    correlationsBarPlot(topPerfs, ds, showPearson, showKendall, showSpearman);
    disp('Correlation plot created.')

    figure;
    rocType = 'standard';
    ROCPlot(topPerfs, ds, rocType, []);
    disp('Standard ROC plot created.')



    figure;
    rocType = 'truncated';
    ROCPlot(topPerfs, ds, rocType, prior);   
    disp('Truncated ROC plot created.')

    figure;
    rocType = 'log-log';
    ROCPlot(topPerfs, ds, rocType, prior);  
    disp('log-log ROC plot created.')

    figure;
    show_llrp = true;
    smoothing = true;
    posteriorAndllrPlusPlot(topPerfs(1), ds, prior, show_llrp, smoothing)
    disp('Posterior and local lr+ plot created.')
   
end
generateTable(Results, tableName);
disp('Table created in ../results')
