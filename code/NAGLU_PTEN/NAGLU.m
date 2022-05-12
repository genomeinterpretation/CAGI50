function out = NAGLU(data_folder)
    %% Generate the reults for NAGLU dataset. 
    % Input
    % data_folder: Path to the data folder.
    % Output
    % out: is a structure with two fields 'perfs' and 'ds'
    % ds: contains the Dataset object for NAGLU, containg dataset specific
    %   information used throughout the pipeline.
    % perfs: a structure array, each entry of which contains the
    %   measures for a method computed using 'perfMetrics()'.
    %   



    data_folder = fullfile(data_folder,'NAGLU_CAGI4');
    data_file = fullfile(data_folder, 'ExperimentalValuesAndPredictionsTop.tsv');
    % Read the files mapping group names to group ids
    group_file = fullfile(data_folder, 'GroupMappings.xls');
    groups = readtable(group_file, 'FileType','spreadsheet', 'ReadVariableNames', true);
    group_ids = cellfun(@(s) s(end-1:end), groups.('Alias'), 'UniformOutput',false);
    group_names = groups.('GroupNameAbbr');

    % Read the file containing the experiemntal values and all predictions.
    data = readtable(data_file, 'FileType', 'text', 'ReadVariableNames',true);
    %disp('NAGLU data size:')
    %disp(size(data))
    fns = data.Properties.VariableNames;
    
    % Extract all experimental values and create all class lables by
    % thresholding the experimental values.
    activity = data.(2);
    allExpVals = activity(~isnan(activity));
    [allClasses, ~] = classfun(activity(~isnan(activity)));
    
    % Create the Dataset object.
    isRegression = true;
    isClassification = true;
    isClinical = true;
    dsDisplayName = 'NAGLU';
    dsAbbr = 'NAGLU';
    expValLabel = 'fractional enzyme activity';
    clsScoreLabel = ['predicted ', expValLabel];
    clsBoundary = 0.15;
    ds = DataSet(allExpVals, allClasses, isRegression, isClassification, isClinical, dsDisplayName, dsAbbr, ...
       expValLabel, clsScoreLabel, clsBoundary);
    out.ds = ds;
    
  
   
    % The predictions are extracted from appropriate columns.
    predictions = data(:,4:length(fns));
    variants = data.(1);
    methods = predictions.Properties.VariableNames;
    perfs={};
    % Iterate over the methods. Map the method to the group name and call
    % the wrapper to evaluate th method.
    for i = 1:length(methods)
        group_id = methods{i}(2:end-2);
        method_id = methods{i}(end);
        ix = find(contains(group_ids, group_id));
        group_name = group_names{ix};
        prediction = predictions.(i);
        ix_nan = isnan(activity) | isnan(prediction);
        activity1 = activity(~ix_nan);
        prediction = prediction(~ix_nan);
        vars = variants(~ix_nan);
        perfMetricWrapper(activity1, prediction, vars, method_id, group_name, false, false)
    end
    
    %--------Experimental Max--------
    activity1 = activity(~isnan(activity));
    stdev = data.(3);
    stdev = stdev(~isnan(activity));
    eMax = ExperimentalMax(activity1, stdev);
    variants = variants(~isnan(activity));
    perfMetricWrapper(activity1, eMax, variants, '' , 'ExpMax', false, true)
    
    %--------Sift and Polyphen2--------
    baseline_file = fullfile(data_folder, 'Baselines.csv');
    baseL_table = readtable(baseline_file, 'FileType', 'text', 'ReadVariableNames',true) ;
    Sift_table =  [baseL_table(: , 'AA_substitution'), baseL_table(: , 'SIFT_score')];
    Sift_table = average_nonunique_variants(Sift_table);
    Pf2_table = [baseL_table(: , 'AA_substitution'), baseL_table(: , 'Polyphen2_HVAR_score')];
    Pf2_table = average_nonunique_variants(Pf2_table);
    Sift_table = innerjoin(data,  Sift_table ,'LeftKeys', 1, 'RightKeys', 1);
    Pf2_table = innerjoin(data,   Pf2_table ,'LeftKeys', 1, 'RightKeys', 1);
    Sift_table = Sift_table(~(isnan(Sift_table.SIFT_score)|isnan(Sift_table.relative_activity)),:);
    Pf2_table = Pf2_table(~(isnan(Pf2_table.Polyphen2_HVAR_score)|isnan(Pf2_table.relative_activity)),:);
    perfMetricWrapper(Sift_table.relative_activity, Sift_table.SIFT_score, Sift_table.AA_substitution, '1', 'Sift', true, false);
    perfMetricWrapper( Pf2_table.relative_activity, 1-Pf2_table.Polyphen2_HVAR_score,  Pf2_table.AA_substitution, '2', 'PolyPhen', true, false);
    
    out.perfs =perfs;
    
    function [classes, class_index] = classfun(activity)
        % Create class labels by thresholding experiemntal values.
        classes = double(activity<0.15);
        class_index = 1:length(activity);
    end

   
    function perfMetricWrapper(activity, predictions, variants, method_id, group_name, isBaseline, isExpMax)
         % Wrapper function to evaluate the methods.
        [cls, cls_index] = classfun(activity);
        scores = -predictions;
        priors = [0.01, 0.1, mean(allClasses)];
        negated = true;
        isMeta = false;
        perf =  PerfMetrics(activity, predictions, cls, scores, negated, cls_index, allClasses, priors, group_name, method_id, ...
            isMeta, isBaseline, isExpMax);
        perfs = [perfs, perf];
    end
    
end


