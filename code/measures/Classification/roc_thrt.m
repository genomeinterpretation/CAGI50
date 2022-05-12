function [auc, fpr, tpr] = roc_thrt(fpr, tpr)
% ROC curve after taking third root of FPR and TPR to enhance the regions
% with small FPR and TPR

tpr = tpr.^(1/3);
fpr = fpr.^(1/3);
auc = trapz(fpr(end:-1:1),tpr(end:-1:1));
end