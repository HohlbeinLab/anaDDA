function [parameters,bootstrapparamstd] = MLEanaDDA(Dlistdata,rangeD,input)
%% Initialize MLE for 

% With input.nofit selected the MLE is skipped and input parameters are
% produced as output

% maxindex = round(max(Dlistdata(2,:).*Dlistdata(1,:))/(rangeD(2)-rangeD(1))*3)+1; % The data size determines the amount
% rangeD = rangeD(1:maxindex);

%% Already determined parameters 
[Dfixed, fitspecies, fixedspecies] = GeneratefixedDdistributions(input.numberofspecies, input.fixedparameters, rangeD,input);


%% MLE
if input.nofit == false
[parameters, ~, bootstrapparamstd] = MLEfitDynamic(Dlistdata,input.numberofspecies,input.fixedparameters,rangeD, Dfixed,fitspecies,fixedspecies,input);
else
parameters = [1- input.fractionB input.koff1_A input.kon1_A input.Dfree_A input.D1_A; input.fractionB input.koff1_B input.kon1_B input.Dfree_B input.D1_B];
parameters = parameters(1:input.numberofspecies,:);
bootstrapparamstd = zeros(size(parameters));
end