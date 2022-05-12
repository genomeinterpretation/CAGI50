function [ix, cnts] = duplicates(x)
% Assumes that x is sorted. Returns the starting index of the duplicate
% groups. And the group sizes.
x = reshape(x, 1, length(x));
dx = x(2:end) - x(1:(end-1));
dup_ix = find(dx == 0);
cnts = [];
ix = [];
if ~isempty(dup_ix)
    ddup_ix = dup_ix(2:end) - dup_ix(1:(end-1));
    breaks = [0, find(ddup_ix>1), length(dup_ix)];
    for i = 1:(length(breaks) - 1)
        b1 = breaks(i)+1;
        b2 = breaks(i+1);
        ix1 = dup_ix(b1);
        ix2 = dup_ix(b2) + 1;
        cnts = [cnts, ix2-ix1+1];
        ix = [ix, ix1];
    end
end
end