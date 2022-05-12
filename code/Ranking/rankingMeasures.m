function rMeasures = rankingMeasures(ds)
%% Returns the measures used for ranking the methods on a dataset
%% based on the applicable analysis types.


    disp([ds.abbreviation, ' is ranked based on:'])
    rMeasures= struct('specs', '', 'ismax', true, 'type', '');
    m = 1;
    if ds.isRegression
        rMeasures(m) = struct('specs', 'kendall', 'ismax', true, 'type', 'Regression');
        rMeasures(m+1) = struct('specs', 'pearson', 'ismax', true, 'type', 'Regression');
        disp('Regression measure: pearson and kendall');
        m = m + 2;
    end
    
    if ds.isClassification
        rMeasures(m) = struct('specs', 'roc auc', 'ismax', true, 'type', 'Classification');
        disp('Classification measure: AUC')
        m = m + 1;
    end
    
    if ds.isClinical
        rMeasures(m) = struct('specs', 'roc auc_trunc', 'ismax', true, 'type', 'Clinical');
        disp('Clinical measure: Truncated AUC');
        m = m + 1;
    end
    disp(' ')
end