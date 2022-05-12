
source_folder = fullfile('../data/AAMWithoutHGMD/');
target_folder = fullfile('../data/AAMWithoutHGMD/');
probable = false;


    
    source_file1 = 'dbnsfp_predictors_on_set1.txt';
    source_file2 = 'dbnsfp_predictors_on_set2.txt';
    target_file1 = 'dbnsfp_predictors_on_set1_strip.txt';
    target_file2 = 'dbnsfp_predictors_on_set2_strip.txt';

    [M1, G1, M2, G2] = groupsAndMethods();

   
    f1_src = fullfile(source_folder, source_file1);
    f1_trgt = fullfile(target_folder, target_file1);
    opts = detectImportOptions(f1_src, 'FileType', 'text', 'Delimiter','\t');
    opts = setvartype(opts, {'chr','pos', M2{:}}, 'char');  
    table1 = readtable(f1_src, opts);
    table1 = table1(:, [{'chr', 'pos', 'ref', 'alt','Ensembl_geneid'}, M2]);
    writetable(table1, f1_trgt,"FileType",'text','Delimiter','\t');
    
    f2_src = fullfile(source_folder, source_file2);
    f2_trgt = fullfile(target_folder, target_file2);
    opts = detectImportOptions(f2_src, 'FileType', 'text', 'Delimiter','\t');
    opts = setvartype(opts, {'chr','pos', M2{:}}, 'char');  
    table2 = readtable(f2_src, opts);
    table2 = table2(:, [{'chr', 'pos', 'ref', 'alt','Ensembl_geneid'}, M2]);
    writetable(table2, f2_trgt,"FileType",'text','Delimiter','\t');
    