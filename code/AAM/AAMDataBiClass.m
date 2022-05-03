function [TT, TTMethods, G] = AAMDataBiClass(folder, withConfidentVariants)
%% The function is similar to AAMData for the most part, except that the returnedonly bi-class
%% variants are contained in the 
[M1, G1, M2, G2] = groupsAndMethods();

if withConfidentVariants
    file1 = 'main_predictors_on_set1.csv';
    file2 = 'dbnsfp_predictors_on_set1.txt';
else
    file1 = 'main_predictors_on_set2.csv';
    file2 = 'dbnsfp_predictors_on_set2.txt';
end
    file1 = fullfile(folder, file1);
    opts = detectImportOptions(file1, 'FileType', 'text');
    opts = setvartype(opts, {'Chrom','Pos'}, 'char');  %or 'char' if you prefer
    table1 = readtable(file1, opts);
    table1 = [table1(:, 1:25), table1(:,32)];
    table1 = table1(filterBasedOnReviewStatus(table1.Review_status),:);
    table1.pph2_prob = averageStrings(table1.pph2_prob);
    [nrow, ncol] = size(table1);

      test_table = unique(table1(:, 1:4));
      [nrow2, ~] = size(test_table);
      if nrow ~= nrow2
          error("need to resolve duplicates in Nilah\'s Table")
      end
    
    file2 = fullfile(folder, file2);
    opts = detectImportOptions(file2, 'FileType', 'text', 'Delimiter', '\t');
    opts = setvartype(opts, {'chr','pos'}, 'char');  %or 'char' if you prefer
    table2 = readtable(file2, opts);
    table2 = table2(:, [{'chr', 'pos', 'ref', 'alt','Ensembl_geneid'}, M2]);
    vars = table2.Properties.VariableNames;
    for ix = 6: length(vars)    
        m = vars{ix};
        if ~isnumeric(table2.(m))
            table2.(m) = averageStrings(table2.(m));
        end
    end
    table2 = average_nonunique_variants(table2, [1,2,3,4],6:length(vars), 5);    
    T = innerjoin(table1, table2, 'LeftKeys', [1,2,3,4],  'RightKeys', [1,2,3,4]);   
    T = T(:,[1:4,27:end]);
    T = outerjoin(table1, T, 'LeftKeys', [1,2,3,4],  'RightKeys', [1,2,3,4]);   
    T = T(cellfun(@(str) ~isempty(str), T.Ensembl_geneid),:);
    %T = T(cellfun(@(str) ~contains(str, remGenes), T.Ensembl_geneid),:);
    G = binGenes(T.Ensembl_geneid);
    %[gene,~,ix_gene] =unique(T.genename);
    ix_path_gene = find(splitapply(@sum, isPathogenic(T.Clinical_significance), G')>=1);
    ix_ben_gene = find(splitapply(@sum, 1-isPathogenic(T.Clinical_significance), G')>=1);
    ix_gene = intersect(ix_ben_gene, ix_path_gene);
    %gene = gene(ix_ben_gene & ix_path_gene);
    ix = arrayfun(@(g) ~isempty(find(ix_gene == g, 1)), G);
    T = T(ix,:);
    TT = T(:, 1:7);
    TT1 = T(:,8:26);
    TT2 = T(:,32:end);
   
    MM1 = TT1.Properties.VariableNames;
    MM2 = TT2.Properties.VariableNames;
    
    for ii = 1:length(M1)
        m = M1{ii};
        if any(strcmpi(m, {'Condel','Turkey','Bologna','pph2_prob'}))
            continue;
        end
        if strcmpi(m, 'PolyPhen')
            v1 = T.PolyPhen;
            v2 = T.Polyphen2_HVAR_score;
        else
            
            match = cellfun(@(m1) firstDisagreement(m,m1), MM1);
            [~,ix1] = max(match);
            mm1 = MM1{ix1};      
            match = cellfun(@(m2) firstDisagreement(m,m2), MM2);
            [~,ix2] = max(match);
            mm2 = MM2{ix2};
            disp([mm1,' ', mm2])
            v1 = T.(mm1);
            v2 = T.(mm2);
        end
       
        disp(m)
        mn = mean(abs(v1-v2), 'omitnan');
        disp(mn)
        disp(std(abs(v1-v2), 'omitnan'))
        p95 = prctile(abs(v1-v2), 95);
        disp(p95);
        if p95 < 0.0001
            
            ix = find(strcmpi(mm2, MM2));
            G2(ix) = [];
            MM2(ix) = [];
        end
    end

    
    TTMethods = [TT1,TT2(:,MM2)];
    G = [G1,G2];
disp('over')
    function i = firstDisagreement(s1, s2)
        i = 0;
        n = min(length(s1), length(s2));
        while (i < n) && (s1(i+1) == s2(i+1))
            i = i+1;
        end
    end
end