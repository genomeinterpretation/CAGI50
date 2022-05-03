function generateTable(Results, fileName)
    %% Call this function to generate an Excel file containing results for a 
    %% group of datasets contained in Results. The excel sheet is stored in the 
    %% results folder of the parent directory. Each sheet in the excel file corresponds
    %% to one dataset.
    % Inputs
    % Results: a struct with the dataset names as fields. Each such field is
    %   itself a struct with two entries 'perfs' and 'ds'. 'perfs' is an
    %   array, each entry of which contains a method's evaluation generated by 
    %   calling 'perfMetics()'. 'ds' is also a struct with fields 'allExpVals' 
    %   (for datasets with regression analysis) and 'allClasses' (for datasets 
    %   with classification analysis). 'allExpVals' and 'allClasses'
    %   contain the all experimental values, after removing nans, and all class 
    %   labels in the dataset. The perfMetrics object might contain fewer 
    %   experimental values and class labels in comparison since a given method
    %   might not make predictions on some variants. 
    % fileName: (string) name of the file where the results should be stored. 
    

    DS = fieldnames(Results);
    % Iterate over the datasets
    for d = 1:length(DS)
        dsAbbr = DS{d};
        %Extract all priors used in the evaluation.
        perfs = Results.(dsAbbr).perfs;
        ds = Results.(dsAbbr).ds;
        allcls = ds.allClasses;
        priors = perfs(1).priors;
        nPriors = length(priors);

        % Create the 'reportedMeasures' structure array containig fields
        % that allow extracting the measures to be reported in the table
        % from 'perfs'. It also contains fields that allows selecting the
        % measures to be displayed, perform statisitcal tests, and
        % appropriate labels to be displayed in the table

        clinicalSpecs = {'clinical tpr', 'clinical fpr', 'clinical LRPlus', 'clinical LRMinus', ...
            'clinical dor', 'clinical mcc', 'clinical ppv',  'clinical ppp'};
        nClinical = length(clinicalSpecs);
        evidences = {'sup', 'mod', 'st'};
        % Entries of 'Specs', 'Priors' and 'Evidence' are used to extract
        % the measures using 'specs2property()' method of 'perfMetrics' class
        Specs = [{'coverage','rSquared', 'rmse', 'pearson', 'spearman', 'kendall', 'roc auc', 'roc auc_trunc' 'roc auc_log'}, ...
            repmat(repmat(clinicalSpecs, 1, nPriors),1, length(evidences))];
        Priors = arrayfun(@(prior) repmat({prior}, 1,length(clinicalSpecs)),priors, 'UniformOutput', false);
        Priors = [Priors{:}];
        Priors = [repmat({nan}, 1, 9), repmat(Priors,1, length(evidences))]; 
        Evidence = arrayfun(@(evi) repmat(evi, 1,length(clinicalSpecs)*nPriors), evidences, 'UniformOutput', false);
        Evidence = [repmat({''}, 1, 9), Evidence{:}];
        
        % 'forRegression', 'forClassification', 'forClinical' allow
        % filtering the measures applicable to a given dataset.
        forRegression = repmat({false}, 1, length(Specs));
        forRegression(2:6) = {true};
        forClassification = repmat({false}, 1, length(Specs));
        forClassification(7) = {true};
        forClinical = repmat({true}, 1, length(Specs));
        forClinical(1:7) = {false};
       
        % 'labelsCol1', 'labelsCol2', 'labelsCol3' contain the measure
        % descriptors displayed in the table.
        labelsCol1 = [repmat({''}, 1, 9), [{'Supporting'},repmat({''}, 1, nClinical*nPriors -1)], ...
            [{'Moderate'}, repmat({''}, 1,nClinical*nPriors -1)], [{'Strong'},repmat({''}, 1,nClinical*nPriors -1)]];
        labelsCol2 = arrayfun(@(prior) [{num2str(prior, 3)}, repmat({''},1,length(clinicalSpecs) -1)], priors, 'UniformOutput', false);
        labelsCol2 = [repmat({''}, 1, 9), repmat([labelsCol2{:}], 1, 3)];
        labelsCol3 = [{'COV', 'R-Squared', 'RMSE', 'Pearson''s correlation','Spearman''s rank correlation', 'Kendall''s Tau', 'AUC', ...
            'Truncated AUC', 'log-log AUC'}, repmat({'TPR', 'FPR', 'LR+', 'LR-', 'DOR', 'MCC', 'PPV', 'PPP'}, 1, nPriors*3)];
        
        % 'StatSig_fcn' contains functions used to compute Statistical
        % significance for some measures.
        rand_auc_log = randomAUC_log(sum(allcls), sum(1-allcls));
        StatSig_fcn = [{''}, {@(bt) bt.p5<0, ''}, {@(bt) bt.p5<0, @(bt) bt.p5<0, @(bt) bt.p5<0, @(bt) bt.p5<0.5, @(bt) bt.p5<0.1, ...
            @(bt) bt.p5<rand_auc_log}, repmat({''}, 1, nPriors*length(clinicalSpecs)*3)];
        
        reportedMeasures = struct('specs', Specs, 'evidence', Evidence, 'prior', Priors, 'labels1', labelsCol1, 'labels2', labelsCol2, ...
            'labels3', labelsCol3, 'statSig_fcn', StatSig_fcn, 'forRegression', forRegression, 'forClassification',  forClassification, ...
            'forClinical', forClinical);

       % Arrange the methods in perfs in the correct order to be displayed
       % in the table. One method is selected from each method group. The
       % representative methods are arranged based on their performance.
       % The baselines and Experimental-Max are displayed towards the end.
        
       perfs = selectMethodsAndBaselines(perfs, ds, length(perfs));

       % Create a matlab Table which can be directly written on to an excel
       % sheet.
        methods = {perfs.method};
        groupLabels = {perfs.group};
        % Create the column headings
        if ds.isClinical
            firstFewHeadings = {'Evidence', 'Prior', 'Measures'};
        else
            firstFewHeadings = {'Measures'};
        end

        T = table('Size', [0, 2*length(perfs)+length(firstFewHeadings)], 'VariableTypes', repmat("string", 1, 2*length(perfs)+length(firstFewHeadings)));
        groupLabels = cellfun(@(g,m) methodDisplayName(g,m), groupLabels, methods, 'UniformOutput', false);
        groupLabelsCI = cellfun(@(str) [str,' CI'], groupLabels, 'UniformOutput', false);
        groupVars = reshape([groupLabels; groupLabelsCI], 1, []);
        T.Properties.VariableNames =[firstFewHeadings, groupVars];

        % Iterate over the measures in 'reportedMeasures', decide if it
        % should be displayed for the given dataset. If yes, then initalize the
        % row of the table with the measure and its 95% confidence interval
        % evaluated for each rported method. 
        i=1;
        for j = 1:length(reportedMeasures)
            measure = reportedMeasures(j);
            %disp(measure);
            if (measure.forRegression && ds.isRegression)|| (measure.forClassification && ds.isClassification) ||...
                    (measure.forClinical && ds.isClinical) ...
                    || strcmpi(measure.specs, 'coverage')
                if ds.isClinical
                    row = {measure.labels1, measure.labels2, measure.labels3};
                else
                    row = {measure.labels3};
                end
                
                try
                    means = arrayfun(@(perf) perf.specs2property(measure, 'val'), perfs);
                catch ME
                    disp('error')
                end
                BT = arrayfun(@(perf) perf.specs2property(measure, 'bt'), perfs);
                if isfield(BT(1), 'p5')
                    try
                        p5 = arrayfun(@(bt) bt.p5, BT);
                    catch ME
                        disp('error')
                    end

                    p95 = arrayfun(@(bt) bt.p95, BT);
                    if strcmpi(class(measure.statSig_fcn), 'function_handle')
                        notSig = arrayfun(@(bt) measure.statSig_fcn(bt), BT);
                        star = [' ', '*'];
                        stars = star((~notSig)+1);
                    else
                        stars = repmat(' ', 1, length(perfs));
                    end
                    entriesMean = arrayfun(@(mean, star) [ num2str(mean,'%.3f'), star], means, stars, 'UniformOutput', false);
                    entriesCI = arrayfun(@(pp5, pp95) [ '[' num2str(pp5,'%.3f'), ', ' num2str(pp95,'%.3f'), ']'], p5, p95, 'UniformOutput', false);
                else
                    entriesMean = arrayfun(@(mean) [ num2str(mean,'%.3f')], means, 'UniformOutput', false);
                    entriesCI = repmat({''}, 1, length(perfs));
                end
                entries = reshape([entriesMean; entriesCI], 1, []);
                row = [row, entries];
                T(i, :) = row;
                i= i + 1;
            end
        end
        % Display the dataset sizes and prior at the bottom of the table.
        
            [nrow, ncol] =size(T);
            T = [T; repmat({''}, 9, ncol)];
            if ds.isRegression && ~ds.isClassification
                T(nrow+5,1:2) = {'Data Set', 'size'};
                T(nrow+6,1:2) = {ds.displayName, num2str(ds.regressionSetSize)};
            elseif ds.classificationSetSize == ds.regressionSetSize || ~ds.isRegression
                T(nrow+5,1:3) = {'Data Set', 'size', 'data prior'} ;
                T(nrow+6,1:3) = {ds.displayName, num2str(ds.classificationSetSize), num2str(mean(ds.allClasses),3)};
            else
                T(nrow+5,1:4) = {'Data Set', 'size regression', 'size classification', 'data prior'};
                T(nrow+6,1:4) = {ds.displayName, num2str(ds.regressionSetSize), num2str(ds.classificationSetSize), num2str(mean(ds.allClasses),3)};
            end
            %T(nrow + 8,1) = {TableText(ds)};
       
        writetable(T, fullfile('../results', fileName), 'FileType', 'spreadsheet', 'Sheet', dsAbbr)

    end
end