function [parameters,bootstrapparamstd,locerrorparam] = MLEanaDDA(Dlistdata,input)
%% Initialize MLE for 

% With input.nofit selected the MLE is skipped and input parameters are
% produced as output

%% Already determined parameters 
%[Dfixed, fitspecies, fixedspecies] = GeneratefixedDdistributions(input.numberofspecies, input.fixedparameters, rangeD,input)