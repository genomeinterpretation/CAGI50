function [aullroc, lfpr, ltpr] = roc_log(fpr, tpr, n_neg, n_pos)
%% Compute the log-log ROC curve from the standard FPR, TPR and the number of positves and negatives.
% Inputs:
% fpr, tpr: the standard fpr and tpr vectors isorted in descending order.
% n_neg, n_pos: the number of negatives and positives in the entire
%   classification dataset. When a method makes predicts only on a subset of the
%   entire classification dataset, for a fair comparison with other
%   methods, n_neg and n_pos should be the class counts on the entire
%   classification data, instead of that on the subset where the method
%   makes predictions.



% The small TPR and FPR values need to be handled appropriately when converting
% to log scale, since the log of a small number approaches -inf. To this end, 
% the idea is to restrict the log-log ROC curve to a fixed range on both the log(TPR)
% and the log(FPR) axes, where the range is a function of n_pos and n_neg. Precisely,
% The log(TPR) axis is restricted to [log(1/n_pos),log(1)]. The log(FPR) axis is 
% restricted to [log(1/n_neg) - log(2), log(1)] (see comments below for the
% rational behind the log(2) term in the log(FPR) range). The dependence on 1/n_neg
% and 1/n_pos comes from the fact that, in absence of duplicate scores, they give 
% the minimum non-zero FPR and TPR values achieved on a dataset with n_neg 
% negatives and n_pos positives (note that the FPR = 0 or TPR = 0 can still be observed).
% Also note that when a method does not predict on the entire data, the TPR, FPR 
% computed on the reduced data can have larger non zero minimum FPR and TPR as the number
% of postives and negatives could be smaller. However, for a fair
% comparison between methods the range of the log-log ROC curve should not
% depend on the method (class counts of the subset on which the method
% makes a prediction). Thus n_pos and n_neg should be the class counts on
% the entire classification data and not reduced to the subset where the
% method makes a prediction

fpr_min = 1/n_neg;
tpr_min = 1/n_pos;

ix = fpr >= fpr_min;
tpr(tpr<tpr_min) = tpr_min;
% When TPR > 1/n_pos is achieved at FPR < 1/n_neg, the minimum TPR value
% achieved is set to the maximum TPR acheived at all FPR < 1/n_neg, otherwise 
% it is 1/n_pos. 
tpr_min_achieved = max(tpr(~ix));
if isempty(tpr_min_achieved)
    tpr_min_achieved = tpr_min;
end
log = @log10;


lfpr = log(fpr);
ltpr = log(tpr);
% Minimum FPR on the log scale is set to log(1/n_neg) - log(2).
% log(2) is the distance between the log transformed smallest and the second 
% smallest non-zero FPR value; i.e., log(2) = log(2/n_neg)-log(1/n_neg)
% The log(TPR) value on the log(2) length FPR segment below log(1/n_neg) is set
% to the log of the minimum TPR value achieved. This adjustment serves to account 
% for the standard ROC curve below FPR of 1/n_neg.
lfpr_min = log(fpr_min)-(log(2*fpr_min)-log(fpr_min));
lfpr(~ix) = lfpr_min;
ltpr_min = log(tpr_min_achieved);
ltpr(~ix) = ltpr_min;
% Compute maximum area under the log-log ROC curve. It is the area of
% a rectangle with sides having length equal to the log(FPR) and log(TPR) range 
max_area =  (log(1)-lfpr_min)*(log(1)-log(tpr_min));
aullroc = trapz(lfpr(end:-1:1),ltpr(end:-1:1)-log(tpr_min))/max_area;
end