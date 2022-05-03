%% Run this script to generate results for NAGLU from raw data
%% The results are stored in '../results'.
delete(gca)
addpathscript;

data_folder = fullfile('..','data');
results_folder = fullfile('..', 'results');

Results = struct();

Results.NAGLU = NAGLU(data_folder);
Results.PTEN = PTEN(data_folder);



file_path = fullfile(results_folder , 'NAGLU_PTEN');
save([file_path,'.mat'], 'Results', '-v7.3');
    

