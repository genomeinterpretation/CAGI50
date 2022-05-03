function [auc, fpr, tpr] = roc_trunc(fpr, tpr, fpr_ub)
%% Computer ROC curve trucated below FPR value of fpr_ub from the standard FPR, TPR values.
if nargin < 3
    fpr_ub = 0.2;
end
ix = fpr>fpr_ub;
% Compute the minimum and maximum (FPR, TPR) values achieved in the data above 
% and below FPR = fpr_ub.
tpr_up = min(tpr(ix));
fpr_up = min(fpr(ix));
tpr_low = max(tpr(~ix));
fpr_low = max(fpr(~ix));
if isempty(tpr_low)
    tpr_low = 0;
    fpr_low = 0;
end
% Interpolate the TPR value when the data doesn't give FPR exactly equal to
% fpr_ub
slope = ((tpr_up-tpr_low)/(fpr_up-fpr_low));
tpr_max = tpr_low + slope * (fpr_ub-fpr_low);
tpr(ix) = tpr_max;
fpr(ix) = fpr_ub;
% Compute the area under the curve and normalize it by dividing by the maximum 
% possible area, which is equal to fpr_ub.
auc = trapz(fpr(end:-1:1),tpr(end:-1:1))/fpr_ub;
end

