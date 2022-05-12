function out = AAM(data_folder, withProbableVariants, database, biClass, dsDisplayName, dsAbbr)
    %% Function to evaluate perfermance on all Annotate all Missense datasets.
    % Input
    % data_folder: the path to the data folder.
    % withProbableVariants: (boolean) set to true to read files containing 
    %   variants with probable and confident annotations. Set to false to read files containing 
    %   variants with confident annotations only.
    % database: (string) Takes value one of the following values
    %   1)'All': to use ClinVar and HGMD variants. 
    %   2)'CV': to only use ClinVar variants
    %   3)'HGMD': to obtain pathogenic annotations only from HGMD
    %   Note that since HGMD is a private dataset, the HGMD variants have
    %   been removed from the data files. Consequenlty 'All' and 'HGMD' options
    %   cannot be used.
    % biClass: (boolean) set to true if only evaluation is to be perfomed
    %   on bi-class variants only.
    % dsDisplayName: (string) containing dataset display name.
    % dsAbbr: (string) containing dataset abbreviation.
    disp(['Evaluating ', dsDisplayName]);
   
    data_folder = fullfile(data_folder,'AAMWithoutHGMD/');
    
    % T and TMethods are tables with same number of rows.
    % T contains the class labels and other annotations for the variants.
    % Each column of TMethods contains the predictions from a method
    % G is cell of strings, giving the group names for the methods 
    % (one entry for each column of TMethods)
    
    [T, TMethods, G] = AAMData(data_folder, withProbableVariants, biClass);
   
    if strcmpi(database, 'CV')
        % Get index of all ClinVar variants.
        ix = isClinvar(T.Clinical_significance);
    elseif strcmpi(database, 'HGMD')
        % Get index of all HGMD variants and the benign/likely benign
        % variants from ClinVar
        ix = ~isClinvarPathogenic(T.Clinical_significance);
    else
        [nrow, ~] = size(T);
        ix = 1:nrow;
    end
    % Filter the tables to contain the correct variants based on the function
    % inputs.
    T = T(ix, :);
    %disp(size(T));
    TMethods = TMethods(ix, :);
    
    % Create the Dataset object.
    allClasses = T.Y;
    isRegression = false;
    isClassification = true;
    isClinical = true;
    expValLabel = '';
    clsScoreLabel = 'classification score';
    clsBoundary = [];
    allExpVals = [];
    ds = DataSet(allExpVals, allClasses, isRegression, isClassification, isClinical, dsDisplayName, dsAbbr, ...
       expValLabel, clsScoreLabel, clsBoundary);
    out.ds = ds;
    [~, ncol] = size(T);
    
    % Iterate over the columns of TMethods. Extract the predictions and
    % call the wrapper to evaluate the performance
    methods = TMethods.Properties.VariableNames;
    perfs = {};
    for i = 1:length(methods)
        m = methods{i};
        scores = TMethods.(m);
        classes = T.Y;
        ix = ~isnan(scores);
        scores = scores(ix);
        classes = classes(ix);
        negated = false;
        roc = ROC(classes, scores);
        % Some method scores have smaller values for pathogenic variants.
        % Detect such methods by computing its AUC. If the AUC is less than
        % 0.5 the method score is negated.
        if roc.auc < 0.5
            scores = -scores;
            negated = ~negated;
            disp('AUC expected to be greater than 0.5');
            disp(['negating scores for ',m])
            roc = ROC(classes, scores);
            disp(['AUC after negating ', m, ' scores: ', num2str(roc.auc)])
        end
        %roc = ROC(classes, scores);
        %disp([m, ': ', num2str(roc.auc)])
        isExpMax = false;
        isBaseline = false;
        if contains(G{i}, {'PolyPhen', 'Sift'}, 'IgnoreCase', true)
            isBaseline = true;
        end
        isMeta = false;
        if contains(G{i},{'Revel', 'MetaLR', 'MetaSVM', 'M_CAP'}, 'IgnoreCase', true)
            isMeta = true;
        end
        % Call the wrapper tp evaluate performance.
        perfMetricWrapper(classes, scores, G{i}, m, negated, isMeta, isBaseline, isExpMax);
    end
    
    out.perfs = perfs;


    function perfMetricWrapper(classes, scores, group, method, negated, isMeta, isBaseline, isExpMax)
        % Wrapper function to evaluate performance.
%         disp('AAM positive:');
%         disp(sum(classes));
%         disp('AAM negatives:');
%         disp(sum(1-classes));
        class_index = 1:length(classes);
        priors = [0.01, 0.1];
        perf =  PerfMetrics([], [], classes, scores, negated,  class_index, allClasses, priors, group, method, isMeta, isBaseline, isExpMax);
        perfs = [perfs, perf];
    end

end





