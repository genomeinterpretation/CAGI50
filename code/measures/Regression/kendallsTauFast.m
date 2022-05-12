function [tau, Q, Tx, Ty, Txy] = kendallsTauFast(x,y, opts)
%% Fast computation of Kendall's Tau-b in O(nlog(n)) [1].
% [1] Knight, William R. "A computer method for calculating Kendall's tau with 
% ungrouped data." Journal of the American Statistical Association 61.314 (1966): 436-439.
% 
% Inputs:
% x and y are vecotrs of same length for which Kendall's Tau is computed.
% opts is struct with boolean field 'presorted'. Set to true, if x is sorted
% and y is reordered to match x, otherwise set to false.
%
n = length(x);
if nargin < 3 || ~isfield(opts, 'presorted') || ~opts.presorted
    [x,ix] = sort(x);
     y = y(ix);
end

% Find ties in x and count the sizes of each set of ties. Tx contains the
% counts.
[ix, Tx] = duplicates(x);
% Find ties in (x,y) as a pair and count the sizes of each set of ties.
% Here an index is tied with another, when both x and y values are identical 
% for the two indices. Txy contains the counts.
Txy = [];
for ii = 1:length(ix)
    ix1 = ix(ii);
    ix2 = ix1 + Tx(ii) - 1;
    ysort = sort(y(ix1:ix2));
    y(ix1:ix2) = ysort;
    [~, txy] = duplicates(ysort);
    Txy = [Txy, txy];
end

 % Now sort the elements of y by iterative merge sorting. On the ith iteration 
 % pairs of consequtive windows of size 2^(i-1) are merge sorted. As a result
 % at the end of the ith iteration, windows of size 2^i appear in sorted order.
 % Repeating this would lead to the entire y vector getting sorted.
 % This iterative merge sorting based approach allows counting the number of
 % discordant pairs as the total number of positions (Q) moved by the elements
 % of y to th left during the merge sort. Note that y is supposed to be ordered 
 % as per the sort order of x before starting the merge sort.
 Q = 0;
 w = 2; 
 % Use variable yy to store the intermediate results inside the loop so that 
 % the elements of y ,not yet placed at the correct position, do not get over 
 % written.
 yy = y;
 while w < 2*n
    for j = 1:w:n
        % During this iteration the window of size w starting at index j
        % gets sorted and stored in yy.
        % j1 and u1 are the first and last index of the w/2 length window to the left.
        % u1 is kept constant inside this loop, whearas j1 is updated. Note
        % that the length could be smaller than w/2 for the last window.
        j1 = j;
        u1 = min(j + w/2 -1,n);
        % j2 and u2 are the first and last index of the window to the right.
        % u2 is kept constant inside this loop, whearas j2 is updated.
        % Note that the length could be smaller than w/2 for the last window.
        j2 = u1 + 1;
        u2 = min(j + w - 1, n);
 
        if u1 == n
            % when the right window is empty, the entries j:n in y is
            % already sorted
            yy(j:n) = y(j:n);
            break;
        end
        % The elements of the two consequtive w/2 length windows are merge 
        % sorted in the loop below. Each position of the merged window in yy is
        % updated one after the other from the smallest to the largest.
        for kk = j:u2
            if (j1 > u1)
                % If all entries of the left window have been placed at
                % the correct position, then remaining entries of the right 
                % window are already in the correct position and are copied to
                % the end of the merged window in yy.
                % None of the entries are moved to the left, so Q is not
                % updated.
                yy(kk:u2) = y(j2:u2);
                break;
            end
            if (j2 > u2)
                % If all entries of the right window have been placed at
                % the correct position, then remaining entries of the left 
                % window can be directly placed at the end of the merged window
                % in yy. Here the entries only move the right so Q is not
                % updated.
                yy(kk:u2) = y(j1:u1);
                break;
            end
            if y(j1) <= y(j2)
                % If the current entry of the left window is smaller than
                % the current entry of the right window, then left window
                % entry is placed at the current positon of the merged
                % window in yy. The entry moves to  the right so Q is not
                % updated. Now that j1 index entry is at its correct
                % position in the merged window, j1 is incremented to the
                % point to the next entry.
                yy(kk) = y(j1);
                j1 = j1 + 1;
            else
                % If the current entry of the left window is larger than
                % the current entry of the right window, then right window
                % entry is placed at the current positon of the merged
                % window in yy. The entry moves to the left so Q is
                % incremented by the number of postions moved. 
                % Now that j2 index entry is at its correct
                % position in the merged window, j2 is incremented to the
                % point to the next entry
                yy(kk) = y(j2);
                Q = Q + j2 - kk;
                j2 = j2 + 1;
            end
        end
    end
    % Now that all pairs of windows are merged, update y so that it can be
    % used in the next iteration.
    y = yy;
    % Since windows of size w in y are sorted, double the merged window
    % size for the next iteration.
    w = w * 2;
 end
% Find ties in y and count the sizes of each set of ties. Tx contains the
% counts. 
[~, Ty] = duplicates(y);

% Total number of pairs of points without accounting foe duplicates.
nPair = n*(n-1)/2;
% Total number of pairs having equal x in both points of the pair.
nPair_x = sum(Tx.*(Tx-1))/2;
% Total number of pairs having equal x in both points of the pair.
nPair_y = sum(Ty.*(Ty-1))/2;
% Total number of pairs having equal (x,y) value for both points of the pair.
nPair_xy = sum(Txy.*(Txy-1))/2;
deno1 = sqrt(nPair - nPair_x);
deno2 = sqrt(nPair - nPair_y);
% effective number of pairs after removing those having either x and/or y
% values identical for both points of the pair.
nPair_effective = nPair -nPair_x - nPair_y + nPair_xy;
% The numerator corresponds to the number of concordant - number of
% discordant pairs. Denominator is to normalize w.r.t the total number of
% pairs after adjusting for ties in x and y.
tau = (nPair_effective - 2*Q)/(deno1*deno2);
end