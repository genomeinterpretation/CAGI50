function out = clinicalMetrics(classes, scores, priors, opts)
    %% Compute all the clinical measures
    % Inputs 
    % classes: the vector of class labels taking value 0 (negatives) or 1 (positives)
    % scores: the vector of a method's prediction scores.
    % priors: is a vector of containing all class priors w.r.t. which all clinical 
    %       measures are computed.
    % opts is a struct that contians some optional inputs containing expensive
    %       computations. If the function is called repeatedly, using opt will be
    %       significantly faster.
    % Fields in opts.
    % EvidenceNotRequired: (bool) true, if clinical evidence thresholds are not 
    %       needed, otherwise false.
    % llrp_xxx: (numeric vector of same length as priors) Here xxx stands for 
    %       sup (supporting), mod (moderate), st (strong), vst (very strong). 
    %       llrp_xxx contains local lr+ threshold at xxx evidence level. ith entry gives 
    %       the local lr+ threshold, at ith prior in 'priors'. 
    % c: (numeric vector of same length as priors) contains the class prior 
    %       dependent constant from the Tavtigan framework
    % fails: (numeric vector of same length as priors) contains the number of 
    %       rules failed for each class prior.
    % llrp_xxx, c and fails can be computed as the output of
    %       prior2thresholds(prior, false, false)
    %
    % Output: is gives as the folowing fields of the 'out' struct
    % scores: input 'scores sorted in ascending order. If the input scores
    %   are discrete (only taking values 1 or 0) the output score is the
    %   vector [0,1]
    % classes: the input 'classes' reordered to match the sorted scores
    % priors: the input 'priors' included in the output
    % is_discrete: (bool) true if the scores are discrete, otherwise false
    % tpr, fpr, tnr, fnr, LRPlus, LRMinus, dor, llrPlus, posterior: vectors of the
    %   length length(out.scores), containing class prior independent quanties: 
    %   TPR (true positive rate), FPR (false positive rate), TNR (true negative rate), 
    %   FNR (false negative rate), LR+ (global positve likelihood ratio), 
    %   LR- (global negative likelihood ratio), DOR (dignostic odds ratio), 
    %   lr+ (local likelihood ratio). The field 'posterior' is a class
    %   prior dependent quantity, but only evaluated at the data prior.
    %   Note that this data prior corresponds to the proportion of
    %   positives in 'classes'. When a method makes prediction only on a
    %   subset of the entire classification set, the data prior computed
    %   from 'classes' maybe different than the data prior on the entire
    %   dataset. Thus this 'posterior' cannot be used for a fair comparison 
    %   of methods. The posterior adjusted to the entire classification
    %   data prior is included in the 'posterior_adj' field, provided the
    %   entire data prior is included in 'priors' variable.
    % posterior_adj, RR, ppv, ppp, mcc: cells with one entry per prior in
    %   'priors'. Each entry is a class prior depndent vector quantity of
    %   length, length(out.scores). 'posterior_adj' contains posterior
    %   adjusted to the priors in 'priors'. Similarly, 'RR', 'ppv', 'ppp' and 
    %   'mcc' contains Relative Risk, Precision, Proportion of positive
    %   predictions, and Matthews correlation coeff, respectively, 
    %   adjusted to the priors in 'priors'.
    % thr_xxx, llrp_xxx, post_xxx, ix_xxx: (vector of length equal to length(priors)) 
    %   Here xxx stands for different evidence levels: sup (supporting), 
    %   mod (moderate), st (strong), vst (very strong).
    %   These are the clinical threshods for scores (thr_xxx), local lr+ (llrp_xxx), 
    %   posterior (post_xxx), index of the score threshold in scores (ix_xxx). 
    % tpr_xxx, fpr_xxx, tnr_xxx, fnr_xxx, ppv_xxx, ppp_xxx, mcc_xxx LRPlus_xxx, 
    %   LRMinus_xxx, dor_xxx: (vector of length equal to length(priors))
    %   Quantities evaluated at the clinical thresholds. 
    % c, fails: (vecors of length equal to length(priors)) 'c' contains the class prior 
    %   dependent constant from the Tavtigan framework; 'fails' gives the number of 
    %   ACMG/AMP rules failed for each class prior.

    if nargin < 4
        opts = struct();
    end
    
    INF = 1000;
    if isfield(opts, 'INF')
        INF = opts.INF;
    end

    classes = double(classes);
    
    priors = reshape(priors, 1, length(priors));
    out.priors = priors;
    n = length(classes);
    classes = reshape(classes,1,n);
    scores = reshape(scores,1,n);
    n_pos = sum(classes);
    n_neg = n - n_pos;
    if n_pos < 10
        warning('Unreliable estimates: Number of positives less than 10.')
    end
    if n_neg < 10
        warning('Unreliable estimates: Number of negatives less than 10.')
    end
    

    scores_u = unique(scores);
    if length(scores_u) < 2
        error('scores has only one unique value');
    end
    
    % If 1 and 0 are the only unique values in the scores, then set the
    % is_discrete flag to true, otherwise set to false.
    is_discrete = false;
    if length(scores_u) == 2
        is_discrete =  all(scores_u == 1 | scores_u == 0);
        if ~is_discrete
            warning('scores have two unique values, but they are not 0 or 1. For discrete computation change the scores to 0 and 1.')
        end
    end
    out.is_discrete = is_discrete;
    
    % Sort the scores and reorder classes if they are not already sorted.
    if ~isfield(opts, 'presorted') || ~opts.presorted
        [scores, ix_sort] = randomSort(scores);
        classes = classes(ix_sort);
    end
    out.classes = classes;
    out.scores = scores;
    % Compute TPR and FPR
    out.tpr = double(cumsum(classes, 'reverse'))/n_pos;
    out.fpr = double(cumsum(1-classes, 'reverse'))/(n-n_pos);
    % For duplicate scores, the highest TPR and FPR achieved on the
    % duplicates is used.
    % since the scores are sorted the duplicates appear at consecutive
    % indices. ixx(i) contains the starting index of the ith set of duplicates. 
    % cnt(i) is the number of duplicates in the ith set of duplicates.
    [ixx, cnt] = duplicates(scores);
    for i = 1:length(ixx)
        ix_i = ixx(i);
        cnt_i = cnt(i);
        % the highest fpr and tpr values for the ith set of duplicates are
        % at ix_i
        out.tpr(ix_i:ix_i+(cnt_i-1)) = out.tpr(ix_i);
        out.fpr(ix_i:ix_i+(cnt_i-1)) = out.fpr(ix_i);
    end
    
    % When the scores contain only 0 and 1, the vectors containing quantities evaluated 
    % at each score value are reduced to a 2 dimensional vector containing the values at 
    % score 0 and score 1. This is done for TPR, FPR, TNR, FNR, llrPlus, posterior, 
    % DOR, LRPlus, LRMinus, PPP, PPV, MCC. The scores returned with the out struct is also 
    % converted to the vector [0, 1]
    if is_discrete
        out.scores = [0,1];
        n=2;
    end
    % 
    if is_discrete
        tpr0 = out.tpr(find(scores==0,1));
        tpr1 = out.tpr(find(scores==1,1));
        fpr0 = out.fpr(find(scores==0,1));
        fpr1 = out.fpr(find(scores==1,1));
        out.tpr = [tpr0, tpr1];
        out.fpr = [fpr0, fpr1];
    end
    
    % Compute TNR, FNR, DOR, LRPlus and LRMinus from TPR and FPR. 
    % Handle infinte values and 0/0 when they occur.
    out.tnr = 1 - out.fpr;
    out.fnr = 1 - out.tpr;

    out.LRPlus = out.tpr./out.fpr;
    out.LRMinus = out.fnr./out.tnr;
    out.LRMinus(isnan(out.LRMinus)) = out.LRMinus(find(~isnan(out.LRMinus), 1, 'first')); 
    %out.dor = (out.tpr.*out.tnr)./(out.fpr.*out.fnr);
    out.dor = out.LRPlus./out.LRMinus;
    out.LRPlus = makefinite(out.LRPlus);
    out.LRMinus = makefinite(out.LRMinus);
    out.dor = makefinite(out.dor);

    
    % Compute the posterior assuming the class prior is same as the data
    % prior.
    if ~is_discrete
        % Compute the posterior by binning the scores.
        posterior = posteriorWindow(classes, scores);
    else
        % Compute the posterior as the class proportions at scores 0 and 1.
        posterior = nan(1, 2);
        posterior(1) = mean(classes(scores==0));
        posterior(2) = mean(classes(scores==1));
    end
    out.posterior = posterior;
    data_prior = mean(classes);
    % Convert the posterior to local lr+. Handle infinte values and 0/0.
    out.llrPlus = makefinite((posterior./(1-posterior))/(data_prior/(1-data_prior)));
    out.RR = arrayfun(@(prior) locallrPlus2RR(out.llrPlus, prior), priors, 'UniformOutput', false);
    % Convert local lr+ to the posteriors adjusted to the class priors
    % in 'priors'.
    out.posterior_adj = arrayfun(@(prior) locallrPlus2Posterior(out.llrPlus, prior), priors, 'UniformOutput', false);
    
    % Compute the class prior dependent measures (PPP, PPV, MCC) at all priors in 'priors'
    [out.ppp,out.ppv, out.mcc] = arrayfun(@(prior) priorDependentMetric(prior), priors, 'UniformOutput', false);
    
    priors = reshape(priors, 1, length(priors));
    
    % Computing the clinical thresholds.
    if ~isfield(opts, 'EvidenceNotRequired')|| ~opts.EvidenceNotRequired
        % Computing the local lr+ thresholds is expensive. When the function is called repeatedly, 
        % it is advised to give the thresholds as input to the function as fields in opt.
        % If the local llr+ thresholds are not given in opt, compute it.
        if ~isfield(opts, 'llrp_sup')
            [out.llrp_sup, out.llrp_mod, out.llrp_st, out.llrp_vst, out.c, out.fails] = ...
                arrayfun(@(prior) prior2lrThresholds(prior, false, false), priors);
        else
           out.llrp_sup = opts.llrp_sup; out.llrp_mod = opts.llrp_mod; out.llrp_st = opts.llrp_st; out.llrp_vst = opts.llrp_st;
           out.c = opts.c;
           out.fails = opts.fails;
        end
        % Convert the local lr+ thresholds to class prior adjusted
        % posterior threshods.
        out.post_sup = arrayfun(@(t, prior) locallrPlus2Posterior(t, prior), out.llrp_sup, priors);
        out.post_mod = arrayfun(@(t, prior) locallrPlus2Posterior(t, prior), out.llrp_mod, priors);
        out.post_st = arrayfun(@(t, prior)  locallrPlus2Posterior(t, prior), out.llrp_st, priors);
        out.post_vst = arrayfun(@(t, prior) locallrPlus2Posterior(t, prior), out.llrp_vst, priors);
    
        % Derive the clincal score thresholds from the local lr+ thresholds
        % Also compute TPR, FPR, PPV, PPP, MCC, global LR+, lR- and DOR at
        % the score thresholds.
        [out.fpr_sup, out.tpr_sup, out.ppv_sup, out.ppp_sup, out.LRPlus_sup, out.LRMinus_sup, out.dor_sup, out.mcc_sup, ...
            out.ix_sup, out.thr_sup] = arrayfun(@(prior, t) scoreThresholds(out.llrPlus, prior, t), priors, out.llrp_sup);
        [out.fpr_mod, out.tpr_mod, out.ppv_mod, out.ppp_mod, out.LRPlus_mod, out.LRMinus_mod, out.dor_mod, out.mcc_mod,...
            out.ix_mod, out.thr_mod ] = arrayfun(@(prior, t) scoreThresholds(out.llrPlus, prior, t), priors, out.llrp_mod);
        [out.fpr_st, out.tpr_st, out.ppv_st, out.ppp_st, out.LRPlus_st, out.LRMinus_st, out.dor_st, out.mcc_st,...
            out.ix_st, out.thr_st] = arrayfun(@(prior, t) scoreThresholds(out.llrPlus, prior, t), priors, out.llrp_st);
        [out.fpr_vst, out.tpr_vst, out.ppv_vst,out.ppp_vst, out.LRPlus_vst, out.LRMinus_vst, out.dor_vst, out.mcc_vst, ...
            out.ix_vst, out.thr_vst] = arrayfun(@(prior, t) scoreThresholds(out.llrPlus, prior, t), priors, out.llrp_vst);
    
    end

    if (isinf(out.dor_mod))
        disp('here')
    end



    function [fpr, tpr, ppv, ppp, LRPlus, LRMinus, dor, mcc, ix_thr, thr] = scoreThresholds(llrp, prior, t)
        % get index of all scores with local lr+ greater than the lr+
        % threshold
        ix = find(llrp>=t);
        % Compute the fraction of scores above an index in ix that have
        % local lr+ above the threshold. And find that smallest index (ix_thr) such
        % that all scores above it have local lr+ above the threshold.
        confvec = (length(ix):-1:1)./(n-ix+1);
        ix_conf = confvec == 1.0;
        ix_thr = min(ix(find(ix_conf)));
        
        
        if ~isempty(ix_thr)
            % Extract TPR, FPR, TNR and FNR at ix_thr.
            tpr = out.tpr(ix_thr);
            fpr = out.fpr(ix_thr);
            tnr = out.tnr(ix_thr);
            fnr = out.fnr(ix_thr);
            
            % Extract PPP, PPV, MCC at ix_thr. Note that these are class prior
            % dependent quantities.
            prior_ix = find(prior==priors);
            pppvec = out.ppp{prior_ix};
            ppvvec = out.ppv{prior_ix};
            mccvec = out.mcc{prior_ix};
            ppp = pppvec(ix_thr);
            mcc = mccvec(ix_thr);
            if tpr == 0
                ppv = 0;
            else
                ppv = ppvvec(ix_thr);
            end

            % Extract global LR+, LR- and DOR at ix_thr.
            if ppp == 0
                LRPlus = 0;
            else
                LRPlus = out.LRPlus(ix_thr);
            end
            if tnr == 0 && fnr == 0
                LRMinus = 0;
            else
                LRMinus =  out.LRMinus(ix_thr);
            end
            if LRPlus == 0 && LRMinus ==0
                dor = 0;
            else
                dor = out.dor(ix_thr);
            end
            
            % Set the clinical threshold to the score achieved at ix_thr.
            thr = out.scores(ix_thr);
            % If the local lr+ at ix_thr is higher than the lr+ threshold,
            % then adjust the score threshold to an interploated value between 
            % out.scores(ix_thr-1) and out.scores(ix_thr). Only do
            % interpolation when scores are continuous.
            if ~is_discrete && ix_thr > 1
                ub = llrp(ix_thr);
                lb = llrp(ix_thr-1);
                slope = (ub-lb)/(thr - out.scores(ix_thr-1));
                if slope > 0 && isfinite(slope)
                    delta = (ub-t)/slope;
                    thr = thr-delta;
                end
            end
        else
            % If local lr+ threshold is not achieved, set the measures to
            % approapriate values.
            tpr = 0.0;
            fpr = 0.0;
            ppp = 0.0;
            LRPlus = 0.0;
            ix_thr = nan;
            ppv = 0.0;
            thr = nan;
            LRMinus = 1;
            dor = 0;
            mcc = 0;
        end
    end

     function x = makefinite(x)
        % Handle infinite values by replacing with maximum or minimum
        % finite values
        maxx = max(x(isfinite(x)));
        minx = min(x(isfinite(x)));
        x(x == inf) = max([maxx, INF]);
        x(x == -inf) = min([minx, -INF]);
    end

    function [ppp, ppv, mcc] = priorDependentMetric(alpha)
        % Compute PPP, PPV, MCC
        ppp = alpha*out.tpr+ (1-alpha)*out.fpr;
        ppv = alpha*out.tpr./ppp;
        ppv = ppv;
        %The formula below is from Ramola et al. Pacific Symposium  on
        %Biocomputing, 2019
        mcc = sqrt((alpha*(1-alpha))./(ppp.*(1-ppp))).*(out.tpr-out.fpr);
        mcc = mcc;
    end

end