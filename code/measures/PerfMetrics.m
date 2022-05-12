classdef PerfMetrics

    % An object of this class stores all the performance measures computed
    % for a method. It also records the true experimental values and their predictions
    % for the regression analysis and the ground truth class labels and the
    % corresponding scores predicted by the method. These are the main 
    % inputs to the class constructor. The inputs for the clinical analysis 
    % is same as that for the classification analysis. For a detailed
    % description of the object's properties see below.
    
    properties
        % Default properties
        nboot = 1000    % number of bootstrap samples. 

        %Properties given as input to the constructor
        expVals         % (numeric n_ev x 1 or empty vector) Contains the n_ev true experimental values
                        % of the varaints, where the method makes a prediction.
                        % Serves as the ground truth for regression type
                        % analysis. Set to an empty vector, if regression analysis is
                        % not applicable.

        predictions     % (numeric n_ev x k or empty vector) Contains the n_ev predicted experimental values 
                        % for the regression analysis. 
                        % Typically, k=1. For Experimental-Max k>1. 
                        % If a regression analysis is not applicable,
                        % set to empty vector.

        classes         % (numeric n_c x 1 or empty vector) Contains the n_c true class labels
                        % for the varaints included in the classification
                        % analysis and where the method makes a prediction.
                        % 1: positives (pathogenic, loss of function, destabilization)
                        % 0: negatives (benign, wild-type function, stable)
                        % Serves as the ground truth for the classification
                        % and the clinical analysis.
                        % If not directly given by the data provider, it is obtained by 
                        % thresholding expVals using a reasonable threshold(s).
                        % The size of classes could be different from the
                        % size of expVals.
                        % If a classification analysis is not applicable,
                        % set to empty vector.

        scores          % (numeric n_c x k or empty vector) Contains the n_c predicted scores for the 
                        % classification analysis and clinical analysis. 
                        % This might be derived from the prediction using reasonable threshold(s) 
                        % or directly given by the data provider. 
                        % Typically, k=1. For Experimental-Max k>1.
                        % If a classification analysis is not applicable,
                        % set to empty vector

        negated         % (boolean or empty vector) If the scores were created from the predictions,
                        % this field indicates if scores = predictions (false) or 
                        % scores = -predictions (true). Typically, the predictions 
                        % should be negated when small values correspond to
                        % class 1. This field is used while plotting to
                        % decide when to unnegate the scores so that the
                        % plots are displayed w.r.t. the original
                        % predictions
                        % If a classification analysis is not applicable,
                        % set to empty vector.

       class_index      % (numeric n_c x 1 or empty vector) 
                        % index array mapping classes and scores to expVals and predictions.
                        % ith entry contains the index of expVal corresponding to ith label in classes 
                        % Used in scatterplots to color the (prediction, expVal) points based on its class label
                        % If a classification analysis is not applicable,
                        % set to empty vector
                        
                        
        allClasses      % (numeric n_ac x 1 or empty vector)  classes only contains the ground truth labels
                        % for the variants where the method makes a prediction. Thus classes 
                        % does not reflect the entire classification set.
                        % allClasses contains all ground truth labels. It
                        % is used to determine a method's coverage, data
                        % prior and the number of positives and negatives
                        % for the log-log ROC computation. For a fair comparison
                        % between two methods based on log-log ROC, allClasses input 
                        % should be identical.
                        % If a classification analysis is not applicable,
                        % set to empty vector


        priors          % (numeric 1 x n_priors or empty vector) priors (proportion of positives in the 
                        % reference population) to be used for the clinical analysis.
                        % Some standard options for the prior are 0.01: for screening, 
                        % 0.1: for dignostic setting data_prior: proportion of positives in
                        % allClasses.
                        % Set to empty vector if clinical analysis is not
                        % applicable.

        group           % (string) The group name or abbreviation of the method. 
                        % Sorting the methods based on their performance requires this field.
                       

        method          % (string) The name or id of a method. 

        isMeta          % (boolean) true if the method is a meta predictor, otherwise false.
        isBaseline      % (boolean) true if the method is a baseline predictor (other then
                        % Experimental-Max), otherwise false.
        isExpMax        % (boolean) true if the method is  Experimental-Max, otherwise false.
 
        
        variants        % (cell array of strings n_ev x 1) (optional) Contains a unique identifier of each variant
                        % where the method makes a prediction. 
                        % Could be the protein/DNA/RNA change string.
                        % The number of entries are same as that in expVals
                        % default: {}

        positions       % (numeric n_ev x 1) (optional) Contains the position of the variants on the protein.
                        % Only variants where the method makes a prediction are included.
                        % The number of entries are same as that in expVals
                        % default: []
        
        %Properties computed from the constructor inputs
        size_reg        % Size of the regression set (n_ev).
        size_class      % Size of the classification set (n_c).
        data_prior      % proportion of positives in allClasses.
        label           % string concatenating id, group and method separated by hyphens.
        coverage        % size of classes divided by size of allClasses
        
        % The properties below contain the various regression, 
        % classification and clinical measures evaluated for the method.
        % Each measure is contained in a struct. For a scalar measure, the
        % struct contains two fields 'val' and 'bt'. For a vector measure
        % the struct contains only one field, 'val'. val contains the 
        % measure evaluated on the entire regression or classification data. bt is a struct
        % containing a summary of the measure evaluated on 1000 bootsrapped samples generated from 
        % the regression or classification data. For the classification data, bootstrapping 
        % is done on the two classes separately. bt has fields 'mean', 'median', 'std', 'p5', 
        % 'p95' giving the mean, median, standard deviation, 5th and 95th percentile, respectively,
        % of the bootstrapped evaluations. For Experimental-Max, instead of bootstrap, 
        % performance is evaluated on the random predictions given as the column
        % of the predictions/scores matrix. For scalar quantities, val is assigned the mean of 
        % these evaluations and bt contains five summary statistic values
        % of the evaluations. For vector quantities, val contians the the
        % evaluation on the first column of the prediction/scores matrix.

        % The following five properties contain the computations for the five regression measures. 
        % Each reagression measure is a scalar quantity. 
        rSquared        % R-Squared 
        rmse            % Root Mean Squared Error
        pearson         % Pearson's correlation
        spearman        % Spearman's correlation
        kendall         % Kendall's Tau

        % roc is a struct containing the evaluations of the classification
        % measures. It contains six ROC curves and the corresponding AUC values
        % given in the following fields: 
        % 1) standard ROC fields: fpr, tpr, auc
        % 2) Truncated ROC fields: fpr_trunc, tpr_trunc, auc_trunc
        % 3) log-log ROC fields: fpr_log, tpr_log, auc_log
        % 4) square-root ROC fields: fpr_sqr, tpr_sqr, auc_sqr 
        % 5) third-root ROC fields: fpr_thr, tpr_thr, auc_thr 
        % 6) fourth-root ROC fields: fpr_for, tpr_for, auc_for 
        % Fields with fpr and tpr prefix represent vector quantities. 
        % Fields with auc prefix represent scalar quantities.
        roc

        % 'clinical' is a struct containing the evaluations of the clinical
        % measures. It contains the four clinical thresholds (supporting, moderate,
        % strong, very strong) defined on the scores, posterior and the local lr+ space. 
        % The thresholds are computed w.r.t. all priors in the 'priors' property. 
        % 'clinical' also contains the measures evaluated at each clinical threshold.
        % Additionlly, it contains quantities such as posterior and local lr+ 
        % evaluated at each score value. A complete list of quantities in
        % 'clinical' along with the field names is given below. 
        % xxx in the list below stands for sup (supporting), mod (moderate), st (strong), vst
        % (very strong). 
        % Some of the fields below are scalar quantities, evaluated separately w.r.t.
        % each class prior in 'priors'. In spite of being a scalar, to store evaluations 
        % for all class priors, their 'val' and 'bt' field contain vectors with one entry 
        % for each class prior. 
        % Some fields below are vector quantities (length same as scores), 
        % evaluated separately w.r.t. each class prior in 'priors'. Their 'val' field is a cell
        % array whose ith cell contains the vector evaluated w.r.t. the ith prior.
        % 1) clinical threshold fileds (scalar evaluated w.r.t. each prior):  thr_xxx (score threshold), 
        %   llrp_xxx (local lr+ threshold), post_xxx (posterior threshold), 
        %   ix_xxx (index of the score threshold in 'scores'). 
        % 2) Quantities evaluated at clinical thresholds (scalar evaluated w.r.t. each prior): 
        %   tpr_xxx (true positive rate), fpr_xxx (false positive rate), tnr_xxx (true negative rate), 
        %   fnr_xxx (false negative rate) ppv_xxx (precision), ppp_xxx (percent of positive predictions), 
        %   mcc_xxx (Matthews correlation coeff.), LRPlus_xxx (global LR+), LRMinus_xxx (global LR-) 
        %   dor_xxx (Diagnostic odds ratio).
        % 3) vector quantities (either independent of the class prior or evaluated w.r.t. data prior only): 
        %   posterior (posterior reflecting the data prior), llrPlus (local lr+), tpr (true positive rate), 
        %   fpr (false positive rate), tnr (true negative rate), fnr (false negative rate), LRPlus (global LR+), 
        %   LRMinus (global LR-) dor (Diagnostic odds ratio).
        % 4) vector quantities (evaluated w.r.t. all priors): 
        %   posterior_adj (posterior adjusted w.r.t. each prior), ppv (precision adjusted w.r.t. each prior), 
        %   ppp (percent of positive predictions adjusted w.r.t. each prior),
        %   mcc (Matthews correlation coeff. adjusted w.r.t. each prior), 
        %   RR (relative risk adjusted w.r.t. each prior).
        % 5) Other: priors ('val' contains the priors used for clinical evaluation, 
        %   same as the 'priors' property of the class object. 'bt' can be ignored), 
        %   scores ('val' contains sorted 'scores'; vector), classes ('val' contains the classes reordered to
        %   match the sorted scores; vector), c ('val' contains the class prior dependent constant from 
        %   the Tavtigan framework; scalar computed w.r.t. each prior), 
        %   fails ('val' gives the number of ACMG/AMP rules failed for each class prior; 
        %   scalar computed w.r.t. each prior),
        clinical
    end
    
    methods
        function obj = PerfMetrics(expVals, predictions, classes, scores, negated, class_index, ...
                allClasses, priors, group, method, isMeta, isBaseline, isExpMax, variants, positions)
            
            if any(any(isinf(expVals))) || any(any(isnan(expVals)))
                error('nan or inf in experimental values')
            end
            if any(any(isinf(predictions))) || any(any(isnan(predictions)))
                error('nan or inf in predictions')
            end
            if any(any(isinf(classes))) || any(any(isnan(classes)))
                error('nan or inf in classes')
            end
            if any(any(isinf(scores))) || any(any(isnan(scores)))
                error('nan or inf in scores')
            end

            obj.expVals = expVals;
            obj.predictions = predictions;
            obj.classes = classes;
            obj.scores = scores;

            obj.size_reg = length(expVals);
            obj.size_class = length(classes);

            if nargin < 15
                positions = [];
            end
            if nargin < 14
                variants = {};
            end

            obj.data_prior = mean(allClasses);
            
            obj.positions = positions;
            obj.variants = variants;
            obj.group = group;
            obj.method = method;
            obj.isMeta = isMeta;
            obj.isBaseline = isBaseline;
            obj.isExpMax = isExpMax;
            obj.label = [group,'-',method];
            obj.priors = priors;
            obj.allClasses = allClasses;
            obj.coverage = length(classes)/length(allClasses);
            if obj.coverage > 1
                error('Coverage cannot be greater than 1')
            end
            obj.class_index = class_index;
            obj.negated = negated;
            
            % Compute measures for regression analysis.
            if ~isempty(expVals)
                obj.rSquared = bootstrap_reg(@rSquared, expVals, predictions, obj.nboot);
                obj.rmse =  bootstrap_reg(@rmse,expVals, predictions, obj.nboot);
                obj.pearson =  bootstrap_reg(@pearson, expVals, predictions, obj.nboot);
                obj.spearman =  bootstrap_reg(@spearman, expVals, predictions, obj.nboot);
                obj.kendall = bootstrap_reg(@kendallsTauFast, expVals, predictions, obj.nboot);
            end
            
            % Compute measures for classification analysis. To ensure that
            % the log-log ROC among methods with different coverages can be
            % compared fairly, the number of positive and negative
            % examples in the entire data, extracted from allClasses, should be 
            % given to the ROC function.
            if ~isempty(classes)
                optsROC.n1 = sum(allClasses);
                optsROC.n0 = sum(1-allClasses);
                obj.roc = bootstrap_cls(@ROC, classes, scores, optsROC, obj.nboot);
            end
            
            % Compute measures for clinical analysis.
            % Precompute the local lr+ thresholds for each prior using
            % Tavtigan's framework. Also obtain the Tavtigan's constant c
            % and the number of ACMG/AMP rules failed.
            if ~isempty(priors) && ~isempty(classes)
                [opts.llrp_sup, opts.llrp_mod, opts.llrp_st, opts.llrp_vst, opts.c, opts.fails] = ...
                    arrayfun(@(prior) prior2lrThresholds(prior, false, false), obj.priors);
                clinicalMetricsFnc = @(cls, scs, optsC) clinicalMetrics(cls, scs, obj.priors, optsC);
                obj.clinical = bootstrap_cls(clinicalMetricsFnc, classes, scores, opts, obj.nboot);
            end
            
        end
        
        function v = specs2property(obj, metric, fld)
            % This function allows convenient acess to the object's properties
            % programmatically, when the properties are provided as string
            % specifiers.
            % fld: should be a string taking value 'val' or 'bt'. 
            % metric: should be a struct with fields 'specs', 'prior' (optional) and 
            %   'evidence' (optional).
            %   'specs' should be a space delimited string that specfies the
            %    evaluation measure to access. For example for AUC, 
            %    spec = 'roc auc'. 
            % For many clinical measures, the prior and/or the evidence level needs 
            % to be specified. For such measures the 'prior' field of the struct 
            % should contain one of the priors in 'priors' (object property). If the 
            % value of the 'prior' field is not in 'priors' an error message is displayed
            % If the evidence level is also required, 'evidence' 
            % should be set to 'sup' (supporting), 'mod' (moderate), 'st' (strong) or 'vst' 
            % (very strong). 
            % For example to access the clinical score threshold for supporting
            % evidence w.r.t. the diagnostic prior of 0.1, 
            % metric = struct('specs', 'clinical thr', 'prior', 0.1, 'evidence', 'sup')
            fields = split(metric.specs, ' ');
            level1 = fields{1};
            prop = obj.(level1);
            if length(fields)>1
                level2 = fields{2};
                if isfield(metric, 'evidence') && contains(metric.evidence, {'sup', 'mod', 'st', 'vst'})
                       level2 = [level2,'_', metric.evidence];
                end
                if isfield(prop, level2)
                    prop = prop.(level2);
                else
                   error('field does not exist');
                end
            end
            if isfield(prop, fld)
                v = prop.(fld);
                if ~isfield(metric, 'prior') || isnan(metric.prior)
                    return
                end
                ix = find(abs(metric.prior - obj.priors) < 10^-3, 1);
                if isempty(ix)
                    error(['The clinical results were computed w.r.t. to the priors: ', num2str(obj.priors), '. ', ...
                        'None of them mathches the given prior of', num2str(metric.prior), ...
                        '. Please run PerfMetrics with ', num2str(metric.prior), ' as the prior.']);
                end
                if isstruct(v)
                    fields = fieldnames(v);
                    for i = 1:length(fields)
                        field = fields{i};
                        v.(field) = v.(field)(:,ix);
                    end
                else
                    v = v(:,ix);
                end
            else
                v=prop;
            end
        end
    end
end

