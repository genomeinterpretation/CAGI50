function out = PTEN(data_folder)
    %% Generate the reults for PTEN dataset. 
    % Input
    % data_folder: Path to the data folder.
    % Output
    % out: is a structure with two fields 'perfs' and 'ds'
    % ds: contains the Dataset object for NAGLU, containg dataset specific
    %   information used throughout the pipeline.
    % perfs: a structure array, each entry of which contains the
    %   measures for a method computed using 'perfMetrics()'.

    data_folder = fullfile(data_folder,'PTEN_TPMT');
    

    % Read the file containg the experimental values. Remove all TPMT
    % variants. Also remove any synonymous, truncated and wild-type variants. Remove
    % any experiemental value that is nan or negative.
    activity_file = fullfile(data_folder, 'ExperimentalValues.csv');
    activity_table = readtable(activity_file, 'FileType', 'text', 'ReadVariableNames', true);
    activity_table = activity_table(strcmp(activity_table.protein,'PTEN'),:);
    activity_table = activity_table(~strcmpi(activity_table.variant,'WT'),:);
    activity_table = activity_table(cellfun(@(str)str(1)~=str(end), activity_table.variant),:);
    activity_table = activity_table(cellfun(@(str)~strcmpi(str(end),'X'), activity_table.variant),:);
    activity_table = activity_table(activity_table.score>=0,:);
    activity_table = activity_table(~isnan(activity_table.score),:);
    [~,uix] = unique(activity_table.(3));
    activity_table = activity_table(uix, :);
    
    

    % Extract all experimental values and create all class lables by
    % thresholding the experimental values.
    allExpVals = activity_table.score;
    [allClasses, ~] =classfun(activity_table.score, activity_table.score);
    
     % Create the Dataset object
    isRegression = true;
    isClassification = true;
    isClinical = true;
    dsDisplayName = 'PTEN';
    dsAbbr = 'PTEN';
    expValLabel = 'relative protein abundance';
    clsScoreLabel = ['predicted ', expValLabel];
    clsBoundary = [0.4, 1.2];
    ds = DataSet(allExpVals, allClasses, isRegression, isClassification, isClinical, dsDisplayName, dsAbbr, ...
       expValLabel, clsScoreLabel, clsBoundary);
    out.ds = ds;
    
    % Read the files mapping group names to group aliases, giving the name
    % of the folders containing all methods from a groups.
    groups_mapping = readtable(fullfile(data_folder,  'GroupMappings.xls'),'FileType','spreadsheet', 'ReadVariableNames', true);
    group_abbrs = groups_mapping.('GroupNameAbbr');
    group_alias = groups_mapping.('Alias');
    group_ids = extractLastDigits(group_alias);
    
    
    % Get all group folders. Go through each folder iteratively. 
    data_folder = fullfile(data_folder, 'Predictions');
    groups = dir(fullfile(data_folder,'Group*'));
    n_groups = length(groups);
    perfs={};
    for j = 1:n_groups
        group = groups(j);
        group_folder = fullfile(data_folder, group.name);
        % Get all files in the group folder containing a method's
        % predictions
        methods = dir(fullfile(group_folder, '*prediction*'));
        ix = find(strcmpi(group.name, group_alias));
        group_name = group_abbrs{ix};
        %disp(group_name)
        group_id = group_ids{ix};
        
        % Iterate over the methods from a group.
        for i = 1:length(methods)
            method = methods(i);
            method_id = num2str(i);
            prediction_table = readtable(fullfile(group_folder,method.name), 'FileType','text', 'ReadVariableNames', true);
            prediction_table = prediction_table(strcmp(prediction_table.Gene_Symbol,'PTEN'),:);
            [~,uix] = unique(prediction_table.(2));
            prediction_table = prediction_table(uix, :);
            % Join the Experimental value table with the predictor table
            % using the variant sting as the key.
            T = innerjoin(activity_table, prediction_table,'LeftKeys', 3, 'RightKeys', 2);  
            activity = T.score;
            prediction = T.Prediction;
            variant = T.variant;
            positions = T.position;
            % Call the wrapper to evaluate the performance of the method
            perfMetricWrapper(activity, prediction, group_name, method_id, false, false, variant, positions);
        end
    end
    %---------------Experimental-Max-----------
    activity = activity_table.score;
    activity1 = activity(~isnan(activity));
    sd = activity_table.sd(~isnan(activity));
    variants = activity_table.variant(~isnan(activity));
    positions = activity_table.position(~isnan(activity));
    expMax = ExperimentalMax(activity1, sd);
    perfMetricWrapper(activity1, expMax, 'ExpMax', '0', false, true, variants, positions);
    
   
    
%     %--------Polyphen2--------
    baseline_file = fullfile(data_folder, 'PolyPhen2.xlsx');
    baseL_table = readtable(baseline_file, 'FileType', 'spreadsheet', 'ReadVariableNames',true, 'Sheet',1) ;
    baseL_table = baseL_table(strcmpi(baseL_table.Gene, 'PTEN'),:);
    Pf2_table = baseL_table(: , {'variants','pph2_prob'});
    Pf2_table = Pf2_table(~isnan(Pf2_table.(2)),:);
    Pf2_table = average_nonunique_variants(Pf2_table);
    [activity, prediction, variants, positions]=  filter_Baseline_table(Pf2_table);
    perfMetricWrapper(activity, 1-prediction, 'PolyPhen', '0', true, false, variants, positions);
    
 
 out.perfs = perfs;
 
     function [activity, prediction, variant, positions] = filter_Baseline_table(tab)
        T = innerjoin(activity_table, tab,'LeftKeys', 3, 'RightKeys', 1);
        activity = T.score;
        prediction = T.(13);
        variant = T.variant;
        positions = T.position;
    end
     

    function [classes, scores, class_index] = classfun(activity, prediction)
        % Create class labels by thresholding experiemntal values.
        index = 1:length(activity);
        ix_pos = index(activity< 0.4);
        ix_neg = index(activity >=0.4 & activity < 1.2);
        classes = [ones(length(ix_pos),1); zeros(length(ix_neg), 1)];
        class_index = [ix_pos'; ix_neg'];
        scores = -[prediction(ix_pos,:);prediction(ix_neg,:)];
    end



    function perfMetricWrapper(activity, predictions, group,  method, isBaseline, isExpMax, variants, positions)
         % Wrapper function to evaluate the methods
        [classes, scores, class_index]  = classfun(activity, predictions);
%         disp('PTEN positive:')
%         disp(sum(classes))
%         disp('PTEN negatives:')
%         disp(sum(1-classes))
        if length(allClasses)< length(classes)
            error('Coverage cannot be greater than 1')
        end
        priors = [0.01, 0.1, mean(allClasses)];
        negated = true;
        isMeta = false;
        perf =  PerfMetrics(activity, predictions, classes, scores, negated, class_index, allClasses,  priors, group, method, ...
            isMeta, isBaseline, isExpMax, [], positions);
        perfs = [perfs, perf];
    end
end



