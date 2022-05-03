function logAUC = randomAUC_log(n_pos,n_neg)
% Computes the log-log AUC of a random classifier on a dataset with n_pos
% positives and n_neg negatives
fpr_width = log10(2*n_neg);
tpr_width = log10(n_pos);
if fpr_width > tpr_width
    area = 0.5*tpr_width^2;
else
    area = 0.5*fpr_width^2;
    area = area + fpr_width* (tpr_width - fpr_width);
end
logAUC = area/(fpr_width*tpr_width);
end