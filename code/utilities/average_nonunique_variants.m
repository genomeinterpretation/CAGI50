function [table2] = average_nonunique_variants(tab, keyCols, avgCols, otherCols)
[nrow,ncol]=size(tab);
if nargin < 2
    keyCols = 1;
end
if nargin < 4
    otherCols =[];
end
if nargin < 3
    avgCols = (max(keyCols)+1): ncol;
    avgCols = setdiff(avgCols, otherCols);
end


%tab = tab(~isnan(tab.(2)),:);
% if length(keyCols) > 1
%     key = cell(nrow,1);
%     for k = keyCols
%         if isnumeric(tab.(k))
%             key = cellfun(@(k,v) [k,'_',num2str(v)], key, tab.(k), 'UniformOutput',false);
%         else
%             key = cellfun(@(k,v) [k,'_',v], key, tab.(k), 'UniformOutput',false);
%         end
%     end
% else
%     key = tab.(keyCols);
% end 
% [vars, jj, ii] = unique(key);
[~, jj, ii] = unique(tab(:,keyCols));
table2 = tab(jj, :);
[nrow2, ~] = size(table2);
%table2.Properties.VariableNames(1) = tab.Properties.VariableNames(1);
meanFnc = @(x) mean(x, 'omitnan');
for c = avgCols
    fnc = accumarray(ii, tab.(c), [], meanFnc);
    val = fnc(1:nrow2);
    %tt = table(val);
    %table2 = [table2,tt];
    %table2.Properties.VariableNames(c) = tab.Properties.VariableNames(c);
    table2.(c) = val;
end

end
