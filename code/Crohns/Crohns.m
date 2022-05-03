function out = Crohns(data_folder)
  
    data_folder = fullfile(data_folder,'Crohns_CAGI4/');
    
    activity_file = fullfile(data_folder, 'DiseaseLabel.csv');
    activity_table = readtable(activity_file, 'FileType', 'text', 'ReadVariableNames', true);

    allClasses = activity_table.Disease;

     % Create the Dataset object.
    isRegression = false;
    isClassification = true;
    isClinical = false;
    dsDisplayName = 'Crohns';
    dsAbbr = 'Crohns';
    expValLabel = '';
    clsScoreLabel =  'method score';
    clsBoundary = [];
    allExpVals = [];
    ds = DataSet(allExpVals, allClasses, isRegression, isClassification, isClinical, dsDisplayName, dsAbbr, ...
       expValLabel, clsScoreLabel, clsBoundary);
    out.ds = ds;
    
    
    groups_mapping = readtable(fullfile(data_folder,  'GroupMappings.xls'), ...
        'FileType','spreadsheet', 'ReadVariableNames', true);
    group_abbrs = groups_mapping.('GroupNameAbbr');
    group_alias = groups_mapping.('Alias');
    group_ids = extractLastDigits(group_alias);
    data_folder = fullfile(data_folder, 'Predictions');
    groups = dir(fullfile(data_folder, 'Group_*'));
    perfs={};
    for j = 1:length(groups)
        group = groups(j);
        %disp(groups(j));
        group_folder = fullfile(data_folder, group.name);
        methods = dir(fullfile(group_folder, '*prediction*'));
        ix = find(strcmpi(group.name, group_alias));
        abbr = group_abbrs(ix);
        abbr = abbr{1};
        group_id = group_ids{ix};
        %disp(group.name);
        for i = 1:length(methods)
            method = methods(i);
            method_id = num2str(i);
            prediction_table = readtable(fullfile(group_folder,method.name), 'FileType','text', 'ReadVariableNames', true, 'Range', [1,1]);
            prediction_table = prediction_table(~isnan(prediction_table.disease_status),: );
            T = innerjoin(activity_table, prediction_table,'LeftKeys', 1, 'RightKeys', 1);  
            classes = T.Disease;
            scores = T.disease_status;
            %disp([abbr, '-', method_id, ' has scores for ', num2str(length(scores)), ' subjects']);
            isBaseline = false;
            isExpMax = false;
            perfMetricWrapper(classes, scores, abbr, method_id, isBaseline, isExpMax);
        end
    end
   
    out.perfs = perfs;
    
    
    function perfMetricWrapper(classes, scores, group, method, isBaseline, isExpMax)
        disp('Crohns4: positives')
        disp(sum(classes))
        disp('Crohns4: negatives:')
        disp(length(classes)-sum(classes));
        class_index = 1:length(classes);
        priors = 0.013;
        negated = false;
        isMeta = false;
        perf =  PerfMetrics([], [], classes, scores, negated, class_index, allClasses,  priors, group, method, ...
            isMeta, isBaseline, isExpMax, [], []);
        perfs = [perfs, perf];
    end
end





