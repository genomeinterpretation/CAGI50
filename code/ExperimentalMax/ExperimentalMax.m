function emax = ExperimentalMax(y, stdev, K)
% Generates K Experiemntal-Max predictions by sampling from a Gaussian with y
% serving as the mean and stdev as the standard deviation. Here y and stdev 
% correspond to the mean and standard deviation of the replicates. It would 
% be incorrect to use the standard error of the mean, instead of the the 
% standrd deviation of the replicates.
% Input:
% y and stdev: are vectors of same length containg the mean and
%   standard deviation of the replicates in the experiment. y should not
%   contain nan values. Any nan is stdev is interpreted as 0.
% K: is the number of Experimental-Max predictions needed.
% Output:
% emax (numeric vector n x K): each column contains a random Experimental-Max 
%   prediction vector.

    if any(isnan(y))
        error('y contains nans');
    end
    if nargin < 3
        K = 1000;
    end
    n = length(y);
    y = reshape(y, [n,1]);
    stdev = reshape(stdev, [n,1]);
    stdev(isnan(stdev))=0;
    r = randn(n, K);
    emax = y + r.*stdev;
end

