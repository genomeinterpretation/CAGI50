function G = binGenes(genes)

for i = 1 :length(genes)
    g = unique(regexp(upper(genes{i}), ';', 'split'));
    g = cellfun(@(str) str(~isspace(str)), g, 'UniformOutput', false);
    g = g(cellfun(@(str) ~isempty(str), g));
    genes{i} = strjoin(g,';');
end
[g_uni, ~, ix] = unique(genes);
match = diag(ones(length(g_uni),1));
gg = regexp(g_uni, ';', 'split');
for  i = 1:length(g_uni)
    for j = i:length(g_uni)
        gi = gg{i};
        gj = gg{j};
        if any(cellfun(@(g) any(strcmpi(gi, g)),gj))
            match(i,j) = 1;
            match(j,i) = 1;
        end
    end
end

GG = conncomp(graph(match));

G = GG(ix);

end