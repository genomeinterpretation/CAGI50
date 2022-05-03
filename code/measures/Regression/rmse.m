function rmse = rmse(y, prediction)
% Computes the root mean square error of prediction on y.
% prediction and y should be vectors of the same size.
RSS = mean((y-prediction).^2);
rmse = sqrt(RSS);
end
