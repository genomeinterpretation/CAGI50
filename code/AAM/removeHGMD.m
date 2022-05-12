
source_folder = fullfile('../../data/AnnotateEverything/');
target_folder = fullfile('../data/AAMWithoutHGMD/');
probable = false;

if probable
    source_file1 = 'test1.csv';
    source_file2 = 'ANE1.txt';
    target_file1 = 'main_predictors_on_set1.csv';
    target_file2 = 'dbnsfp_predictors_on_set1.txt';
else
    source_file1 = 'test2.csv';
    source_file2 = 'ANE2.txt';
    target_file1 = 'main_predictors_on_set2.csv';
    target_file2 = 'dbnsfp_predictors_on_set2.txt';
end
    f1 = fullfile(source_folder, source_file1);
    table1 = readtable(f1, "FileType",'text');
    opts = detectImportOptions(f1, 'FileType', 'text');
    opts = setvartype(opts, table1.Properties.VariableNames, 'char');  %or 'char' if you prefer
    table1 = readtable(f1, opts);
    table1 = table1(isClinvar(table1.Clinical_significance),:);
    disp('unique clinical significances')
    disp(unique(table1.Clinical_significance))
    disp('main Predictors size')
    disp(size(table1))

    writetable(table1, fullfile(target_folder, target_file1), 'WriteVariableNames', true);
    
    test_table = unique(table1(:, 1:4));
    [nrow , ~] = size(table1);
    [nrow2, ~] = size(test_table);
    if nrow ~= nrow2
        error("need to resolve duplicates in Nilah\'s Table")
    end
    
    f2 = fullfile(source_folder, source_file2);
    table2 = readtable(f2, "FileType",'text', 'Delimiter', '\t');
    opts = detectImportOptions(f2, 'FileType', 'text');
    opts = setvartype(opts, table2.Properties.VariableNames, 'char');  %or 'char' if you prefer
    table2 = readtable(f2, opts);
    [~,ncol2] = size(table2);
    T = innerjoin(table2, table1, 'LeftKeys', [1,2,3,4],  'RightKeys', [1,2,3,4]); 
    disp(T.Properties.VariableNames(ncol2:end))
    disp('unique clinical significances')
    disp(unique(T.Clinical_significance))
    T = T(:, 1:ncol2);
    T.Properties.VariableNames = table2.Properties.VariableNames; 
    disp('dbnsfp predictors size')
    disp(size(T))
    writetable(T, fullfile(target_folder, target_file2), 'WriteVariableNames', true, 'Delimiter', '\t');
    
   