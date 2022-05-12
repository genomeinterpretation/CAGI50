function [auc, fpr, tpr] = roc_sqrt(fpr, tpr)
% ROC curve after taking squareroot of FPR and TPR to enhance the regions
% with small FPR and TPR

tpr = tpr.^(1/2);
fpr = fpr.^(1/2);
auc = trapz(fpr(end:-1:1),tpr(end:-1:1));
end

