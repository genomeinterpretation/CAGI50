%% Generate Annotate all missense results.
%% Since HGMD Data is private, the results are only generated from ClinVar data.
%% The results are stored in '../results'.
% Note that since HGMD variants are not public, the original files containing 
% those variants is not provided in ../data. The files containing only ClinVar 
% variants are provided. 

delete(gca)
addpathscript;
data_folder = fullfile('..', 'data');
results_folder = fullfile('..', 'results');


DispNames = {'Confident variants from ClinVar', 'Confident and probable variants from ClinVar', ...
    'Confident biclass variants from ClinVar', 'Confident and probable biclass variants from ClinVar'};

Results = struct();
% Use only confident variants.
withProbableVariants = false;
biClass = false;
Results.AAM1CV = AAM(data_folder, withProbableVariants, 'CV', biClass, DispNames{1}, 'AAM1CV');
% Use confident and probable variants.
withProbableVariants = true;
Results.AAM2CV = AAM(data_folder,  withProbableVariants, 'CV', biClass, DispNames{2}, 'AAM2CV');
% Use confident biclass gene variants. 
biClass = true;
withProbableVariants = false;
Results.AAM1BiClassCV = AAM(data_folder, withProbableVariants, 'CV', biClass, DispNames{3}, 'AAM1BiClassCV');
% Use confident and probable biclass gene variants.
withProbableVariants = true;
Results.AAM2BiClassCV = AAM(data_folder, withProbableVariants, 'CV', biClass, DispNames{4}, 'AAM2BiClassCV');


file_path = fullfile(results_folder, 'AAM');
save([file_path,'.mat'], 'Results','-v7.3');


