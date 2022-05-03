function [TT, TTMethods, G] = AAMData(folder, withProbableVariants, biClass)
%% Reads the files containing the annotations and predictions.
%% Return two tables 1) TT: has columns identifying the variant, 
%% their class labels and other clinical annotations, and 2) TTMethod:
%% each column of which contains predictions from a method. It also returns 
%% a cell containing the group names of the methods in TTMethods.
% Input:
% folder: the path of the folder containing the AAM files
% withProbableVariants: (boolean) set to true to read files containing 
%   variants with probable and confident annotations. Set to false to read files containing 
%   variants with confident annotations only.
% biClass: (boolean) set to true to return tables with variants from bi-class 
%   genes only, otherwise set to false.
% The predictors are distriubuted in two files (file1 and file2). 
% M1 contains the method names of the first file. 
% M2 contains the method names of the second file.
% G1 and G2 are the group names for the methods in M1 and M2, respectively.
[M1, G1, M2, G2] = groupsAndMethods();

% file1 contains variant descriptors, their class labels, clinical annotations
% and the first set of predictors.
% file2 contains the second set of predictors downloaded from dbnsfp.
% Some predictors in file1 and file2 are same. Only one predictor is 
% kept in that case.
if ~withProbableVariants
    % files containing variants with confident annotations.
    file1 = 'main_predictors_on_set1.csv';
    file2 = 'dbnsfp_predictors_on_set1.txt';
else
    % files containing variants with confident and 
    % probable annotations.
    file1 = 'main_predictors_on_set2.csv';
    file2 = 'dbnsfp_predictors_on_set2.txt';
end
    file1 = fullfile(folder, file1);
    opts = detectImportOptions(file1, 'FileType', 'text');
    opts = setvartype(opts, {'Chrom','Pos'}, 'char');  
    table1 = readtable(file1, opts);
    table1 = [table1(:, 1:25), table1(:,32)];
    % Remove unreliable variant annotations based on review status.
    table1 = table1(filterBasedOnReviewStatus(table1.Review_status),:);
    % PolyPhen-2 predictions appear as string of multiple numbers in the
    % file. Convert these predictions to numeric values by averaging the
    % multiple predictions for a variant.
    table1.pph2_prob = averageStrings(table1.pph2_prob);
    [nrow, ncol] = size(table1);
      test_table = unique(table1(:, 1:4));
      [nrow2, ~] = size(test_table);
      if nrow ~= nrow2
          error("need to resolve duplicates in file 1")
      end
    
    % Read file2 containg the predictors downloaded from dbsnfp
    file2 = fullfile(folder, file2);
    opts = detectImportOptions(file2, 'FileType', 'text', 'Delimiter','\t');
    opts = setvartype(opts, {'chr','pos', M2{:}}, 'char');  
    table2 = readtable(file2, opts);
    table2 = table2(:, [{'chr', 'pos', 'ref', 'alt','Ensembl_geneid'}, M2]);
    vars = table2.Properties.VariableNames;
    for ix = 6: length(vars)    
        m = vars{ix};
        if ~isnumeric(table2.(m))
            % Any predictor with predictions as a string of multiple
            % numbers, is converted to numeric values by averaging the
            % multiple predictions
            table2.(m) = averageStrings(table2.(m));
        end
    end
    % Average the peredictions for rows containg the same variants 
    table2 = average_nonunique_variants(table2, [1,2,3,4],6:length(vars), 5);  
    % Join table1 and table2
    T = innerjoin(table1, table2, 'LeftKeys', [1,2,3,4],  'RightKeys', [1,2,3,4]);   
    T = T(:,[1:4,27:end]);
    % After the first join some variants of table1 not present in table2 would be
    % missing from T. Doing the outer join would add back those variants with table 2 
    % columns having nans for those variants. 
    T = outerjoin(table1, T, 'LeftKeys', [1,2,3,4],  'RightKeys', [1,2,3,4]);   
    
    if biClass
        % Keep variants from non bi-class genes only
        T = T(cellfun(@(str) ~isempty(str), T.Ensembl_geneid),:);
        G = binGenes(T.Ensembl_geneid);
        ix_path_gene = find(splitapply(@sum, isPathogenic(T.Clinical_significance), G')>=1);
        ix_ben_gene = find(splitapply(@sum, 1-isPathogenic(T.Clinical_significance), G')>=1);
        ix_gene = intersect(ix_ben_gene, ix_path_gene);
        ix = arrayfun(@(g) ~isempty(find(ix_gene == g, 1)), G);
        T = T(ix,:);
    end
    % Separate the columns of T 
    % TT contains the varinant descriptors, class labels and other
    % annotations.
    % TT1 contains predicitors from file1 and TT2 contains predictors for
    % file2. 
    TT = T(:, 1:7);
    TT1 = T(:,8:26);
    TT2 = T(:,32:end);
 
    % Detect predictors that are present in both TT1 and TT2.
    % Remove the duplicate predictor from TT2
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
            disp(['Comparing ', mm1,' and ', mm2, ' predictions'])
            v1 = T.(mm1);
            v2 = T.(mm2);
        end
        
        mn = mean(abs(v1-v2), 'omitnan');
        %disp(mn)
        %disp(std(abs(v1-v2), 'omitnan'))
        p95 = prctile(abs(v1-v2), 95);
        %disp(p95);
        if p95 < 0.0001
            ix = find(strcmpi(mm2, MM2));
            disp(['Removing ', mm2, ' from TT2 as it is already contained in TT1'])
            G2(ix) = [];
            MM2(ix) = [];
        end
    end

    
    TTMethods = [TT1,TT2(:,MM2)];
    G = [G1,G2];
%disp('over')
disp('Read data into tables')

    function i = firstDisagreement(s1, s2)
        i = 0;
        n = min(length(s1), length(s2));
        while (i < n) && (s1(i+1) == s2(i+1))
            i = i+1;
        end
    end
end