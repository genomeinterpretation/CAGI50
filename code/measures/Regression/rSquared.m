function rsq = rSquared(y, prediction)
%% Computes the R-squared (coefficeint of determination) to quantify the  
%% performance of prediction on y.
% prediction and y should be vectors of the same size.
TSS = mean((y-mean(y)).^2);
RSS = mean((y-prediction).^2);
rsq = 1 - RSS/TSS;
end

