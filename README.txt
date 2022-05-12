
CAGI Flagship code

Author: Shantanu Jain

The code used for analyzing the datasets in the CAGI Flagship paper is contained in ./code folder. 
Due to data privacy issues only publicly available datasets are provided in this submission.
Based on this criteria, only full datasets for NAGLU and PTEN are provided. The Annotate all
missense data, stripped of all HGMD variants, is provided. Additionally, for the Crohn’s challenge, 
the only class labels and the subject id is provided in the dataset. The data files are contained in ./data folder.

The predictions from only publicly available tools and/or top performing methods in a 
given challenge are provided. The predictions (either merged in the data file or in separate files)
are also contained in ./data folder.

The analysis pipeline can be divided into four parts, 1) Data Interface (DI): reading the experimental 
data and prediction files, 2) Measure Computation (MC): evaluating the predictions by computing relevant
measures, 3) Ranking and selection (RS): ranking methods based on appropriate measures and selecting the 
top performers, along with top performing baseline (if applicable) and Experimental-Max (if applicable)
and 4) Figures and Tables (FT): generating figures and tables. 

The DI code provided in the submission is reduced to only that needed for the provided datasets.  
MC, RS and FT code are provided entirely. 

The perfMetrics class in ./measures is central to the pipeline. An object of this class contains all 
measures (regression and/or classification and/or clinical) computed for a method on a dataset. 
The DI code reads the data and method files and creates a perfMetrics object for each method, where
all relevant measures and their bootstrap summaries are computed with the MC code. The array of 
perfMetrics objects, along with a Dataset object, recording some properties of the dataset, used in 
annotating figures and tables, is stored as a .mat file in the ./results folder. 
The .mat file containing the results is later used by the FT code to generate the figures and tables
after ranking and selecting the methods with RS. The tables are stored in the ./results folder.

In this submission, the scripts for computing the measures and storing them as .mat files are 
1) generateResults_NAGLU_PTEN.m
2) generateResults_AAM.m
3) generateResults_Crohns.m

The scripts to read the generated result files and create the figures and tables is contained in
1) FiguresAndTables_NAGLU_PTEN.m
2) FiguresAndTables_AAM.m
3) FiguresAndTables_Crohns.m

 




