function [auc, fpr, tpr] = roc_fort(fpr, tpr)
% ROC curve after taking fourth root of FPR and TPR to enhance the regions
% with small FPR and TPR

tpr = tpr.^(1/4);
fpr = fpr.^(1/4);
auc = trapz(fpr(end:-1:1),tpr(end:-1:1));
end