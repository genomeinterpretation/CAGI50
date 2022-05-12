%% Run this script to generate results for Crohns from data
%% The results are stored in '../results'.
delete(gca)
addpathscript;

data_folder = fullfile('..','data');
results_folder = fullfile('..', 'results');

Results = struct();

Results.Crohns = Crohns(data_folder);


file_path = fullfile(results_folder , 'Crohns');
save([file_path,'.mat'], 'Results', '-v7.3');
   

