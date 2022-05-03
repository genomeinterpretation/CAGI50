function [perfsByGrp,index_k, index] = firstByGrp(perfs, k)
% pick the first method for eah group. Perfs are assumed to be already sorted
% based on an appropriate ranking measure.
if nargin < 2
    k = length(perfs);
end
groups = {perfs.group};
groups_unique = {};
index = [];
for i = 1:length(groups)
    if ~any(strcmpi(groups_unique, groups(i)))
        index = [index, i];
        groups_unique = [groups_unique, groups(i)];
    end
end
perfsByGrp = perfs(index);
if ~isempty(index)
    perfsByGrp = perfsByGrp(1:min(k,length(perfsByGrp)));
    index_k = index(1:min(k,length(index)));
else
    index_k = [];
end
    
end

