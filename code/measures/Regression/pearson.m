function p = pearson(y, prediction)
%% Computes the pearson correlation coefficient between y and prediction.
% prediction and y should be vectors of the same size.
c = corr([y,prediction], 'type', 'pearson');
p=c(1,2);
end
