function roc = ROC(classes, scores, opts)
%% Computes the six different ROCs and the corresponding AUC values.
% Inputs:
% classes: (numeric vector) contains the class labels taking value 0 (negatives)
%   and 1 (positives).
% scores: contains prediction scores of the methods on average the higher
%   scores should correspond to the class label 1. In some cases the
%   original method predictions needs to be negated for that effect.
% opts: (struct) optional argument having fields 
% 1) 'n0' and 'n1', giving the number of positive and negative points 
%   in the entire classification dat set. These numbers are used to determine 
%   the FPR and TPR range in the log ROC curve. They can be different from 
%   that counted from 'classes' because length of 'classes' could be less than 
%   the size of the entire classification dataset, when a method cannot make
%   predictions on all data points. For a fair comparison between methods
%   based on log-log ROC, the function should be called with the same
%   values of 'n0' and 'n1'. The default values are computed from the
%   'classes'.
% 2) 'fpr_ub': numeric scalar giving the FPR upper bound for Truncated ROC.
%   The default value is 0.2.
% 3) 'preSorted': (bool) This argument can be set to true if the scores are
%   already sorted. This will prevent the function from resorting the
%   scores.
% The six ROCs along with the corresponding variable names are given below. 
% 1) fpr, tpr, auc: The standard ROC
% 2) fpr_trunc, tpr_trunc, auc_trunc: Truncated ROC curve restricted to the 
%   FPR values in the range [0, opts.fpr_ub].
% 3) fpr_log, tpr_log, auc_log: The log-log ROC curve obtained by plotting 
%   log10(TPR) against log10(FPR). See the supplementary material of the 
%   CAGI Flagship paper.
% 4) fpr_sqr, tpr_sqr, auc_sqr: ROC curve obtained by plotting square-root of
%  the standard FPR and TPR. 
% 5) fpr_thr, tpr_thr, auc_thr: ROC curve obtained by plotting third-root of
%  the standard FPR and TPR.
% 6) fpr_for, tpr_for, auc_for: ROC curve obtained by plotting fourth-root of
%  the standard FPR and TPR.

roc = struct();
classes = double(classes);
if nargin < 3 || ~opts.presorted 
    [scores, ix] = randomSort(scores);
    classes = classes(ix);
end
n = length(classes);
n_pos = sum(classes);
n1 = n_pos;
n0 = n-n_pos;
fpr_ub = 0.2;
if nargin >= 3 
    if isfield(opts, 'n1') && isfield(opts, 'n0')
        n1 = opts.n1;
        n0 = opts.n0;
    end
    if isfield(opts, 'fpr_ub') 
        fpr_ub = opts.fpr_ub;
    end
end
classes = reshape(classes,1,n);
roc.tpr = double(cumsum(classes, 'reverse'))/n_pos;
roc.fpr = double(cumsum(1-classes, 'reverse'))/(n-n_pos);
% For duplicate scores, the highest TPR and FPR achieved on the
% duplicates is used.
% since the scores are sorted the duplicates appear at consecutive
% indices. ix(i) contains the starting index of the ith set of duplicates.
% cnt(i) is the number of duplicates in the ith set of duplicates.
[ix,cnt] = duplicates(scores);
for i = 1:length(ix)
    ix_i = ix(i);
    cnt_i = cnt(i);
    % the highest fpr and tpr values for the ith set of duplicates are
    % at ix_i
    roc.tpr(ix_i:ix_i+(cnt_i-1)) = roc.tpr(ix_i);
    roc.fpr(ix_i:ix_i+(cnt_i-1)) = roc.fpr(ix_i);
end
roc.tpr = [roc.tpr,0];
roc.fpr = [roc.fpr,0];
roc.auc =trapz(roc.fpr(end:-1:1), roc.tpr(end:-1:1));
% Compute the log-log ROC and AUC.
[roc.auc_log, roc.fpr_log, roc.tpr_log] = roc_log(roc.fpr, roc.tpr, n0, n1);
% Compute the square-root ROC and AUC.
[roc.auc_sqr, roc.fpr_sqr, roc.tpr_sqr] = roc_sqrt(roc.fpr, roc.tpr);
% Compute the third-root ROC and AUC.
[roc.auc_thr, roc.fpr_thr, roc.tpr_thr] = roc_thrt(roc.fpr, roc.tpr);
% Compute the fourth-root ROC and AUC.
[roc.auc_for, roc.fpr_for, roc.tpr_for] = roc_fort(roc.fpr, roc.tpr);
% Compute the Truncated ROC and AUC.
[roc.auc_trunc, roc.fpr_trunc, roc.tpr_trunc] = roc_trunc(roc.fpr, roc.tpr, fpr_ub);
end