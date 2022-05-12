function [posterior, L ,H] = posteriorWindow(classes, scores, eps2)
% Gives the porortion of positives in a bin containing the score as the posterior.
% Assumes that the scores are sorted in ascending order and the classes are
% ordered accordingly.
classes = double(classes);
scores_pos = scores(classes == 1);
scores_neg = scores(classes == 0);
n = length(classes);
% The range of scores is calculated as the difference between the 95th and
% 5th percentile of the scores. 
range = prctile(scores_pos, 95)-prctile(scores_neg, 5);
% eps is used as half of the bin width.
lambda = 0.05;
eps = lambda*range;
if nargin<3
    % eps2 is used as the minimum number of points around a score to be included 
    % while calculating the posterior at that score.
    % If the number of points in the bin containing a score value is less
    % than eps2 additional points outside the bin are used.
    eps2 =  min(50, max(ceil(2*lambda*n),1));
end

% Handle duplicate scores
newClasses = classes;
[ix, cnt] = duplicates(scores);
for ii = 1:length(ix)
    ix_i = ix(ii);
    cnt_i = cnt(ii);
    newClasses(ix_i: (ix_i+cnt_i - 1)) = mean( classes(ix_i: (ix_i+cnt_i - 1)));
end

low = 1;
high = 1;
posterior = nan(1,n);
L = nan(1,n);
H = nan(1,n);
% Iterate through the scores and compute the posterior at each value.
for i = 1:length(scores)
    % Determine the smallest (l) and the largest (h) score index in the bin around the 
    % current score
    s = scores(i);
    while low<n && scores(low)<s-eps
        low = low + 1;
    end
    while high<=n && scores(high)<=s + eps
        high = high + 1;
    end
    l = low;
    h = high;
    
    % If the number of points in the window is less than eps2 include more
    % points outside the window until eps2 number of points are included.
    while h - l < eps2 && ((h <= n) || (l >1))
        % decrement l if the distance between the lower index score and the
        % current score is smaller than that between the higher index score and the
        % current score, provided the index does not go out of the range,
        % other wise increment h
        if l > 1
            delta_l = abs(scores(l-1)-s);
            if h <= n
                delta_h = abs(scores(h) -s);
                if delta_l < delta_h
                    l = l-1;
                else
                    h = h+1;
                end
            else
                l = l-1;
            end
        else
            h = h+1;
        end
    end
    L(i) = l;
    H(i) = h-1;
    % Calculate the posterior at the ith score as the proportion of positives in the
    % neighborhood around the score.
    posterior(i) = mean(newClasses(l:(h-1)));
end
end

