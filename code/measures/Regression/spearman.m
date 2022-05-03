function p = spearman(y, prediction)
%% Computes the Spearman's correlation between prediction and y.
% prediction and y should be vectors of the same size.
c = corr([y,prediction], 'type', 'spearman');
p=c(1,2);
end