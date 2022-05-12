classdef DataSet
    %% This class contains properties of the dataset


    properties
        allExpVals              % (numeric vector) contains all experimental values excluding nans
        allClasses              % (numeric vector) containing all class labels.
        regressionSetSize       % (numeric) size of the regression set (if applicable), otherwise 0.
        classificationSetSize   % (numeric) size of the classification set, (if applicable), otherwise 0.
        isRegression            % (boolean) true if regression analysis is applicable, otherwise false.
        isClassification        % (boolean) true if classification analysis is applicable, otherwise false.
        isClinical              % (boolean) true if clinical analysis is applicable, otherwise false.
        displayName             % (string) Display name of the dataset to be used in figures and tables.
        abbreviation            % (string) abbreviated name of the dataset to be used as sheet names.
        expValLabel             % (string) The label of the experimental values to be used in the figures, if regression 
                                %   analysis is applicable, otherwise it is ''
        classScoreLabel         % (string) The label of the classification score, if classification analysis 
                                % is applicabe, otherwise it is ''
        classBoundary           % (numeric vector) Contains the class boundary applied to the Experiemntal values
                                %   to create the ground truth labels. In some cases, two boundaries are needed 
                                %   to exclude points not used in classification analysis. In such cases classBoundary
                                %   is of length 2. If the class labels are not defined with such boundaries or
                                %   classification analysis is not applicable, it is an empty vector
                                 
                          
                                
    end

    methods
        function obj = DataSet(allExpVals, allClasses, isRegression, isClassification, isClinical, ...
                displayName, abbreviation, expValLabel, classScoreLabel, classBoundary)
            
            obj.allExpVals = allExpVals;
            obj.allClasses = allClasses;
            obj.regressionSetSize = length(allExpVals);
            obj.classificationSetSize = length(allClasses);
            obj.isRegression = isRegression;
            obj.isClassification = isClassification;
            obj.isClinical = isClinical;
            obj.displayName = displayName;
            obj.abbreviation = abbreviation;
            obj.expValLabel = expValLabel;
            obj.classScoreLabel = classScoreLabel;
            obj.classBoundary = classBoundary;

        end

    end
end