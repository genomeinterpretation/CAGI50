function [x, ix] = randomSort(x)
[nrow, ncol] = size(x);
if ncol ~=1 && nrow ~=1
    error('x should be one dimensional')
end

[x,ix] = sort(x);
[x,ixx] = permuteDuplicatesInSorted(x);
ix = ix(ixx);

% x_u = unique(x);
% dist = pdist2(x_u(:), x_u(:));
% if sum(dist>0) > 0
%     min_d = min(min(dist(dist>0)));
% else
%     min_d = 10^-6;
% end
% nn = length(x);
% %for ii = 1:10
% delta  = rand(nrow, ncol).*(2*double(rand(nrow,ncol)>=0.5) - 1)*10^-1*(min_d/2);
% [~,ix] =  sort(x+delta);
% x = x(ix);
%end
end
