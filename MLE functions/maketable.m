function [tableout,totalparam] = maketable(fittingparameters,fixedparameters,indexfittingparameters,numberofspecies)
totalparam = fixedparameters;
totalparam(indexfittingparameters) = fittingparameters;
totalparam = totalparam(1:numberofspecies,:);
species = 1:length(totalparam(:,1));
species = species';
fraction = totalparam(:,1);
fraction(1) = 1-sum(fraction(2:end));
koff = totalparam(:,2);
kon = totalparam(:,3).*koff;
Dfree = totalparam(:,4);
if sum(totalparam(:,5))>0
    D1 = min(totalparam(:,4),totalparam(:,5));
    D2 = max(totalparam(:,4),totalparam(:,5));
tableout = table(species,fraction,koff,kon,D1,D2);
    totalparam(:,5) = D1;
    totalparam(:,4) = D2;
    totalparam = totalparam(:);
else
tableout = table(species,fraction,koff,kon,Dfree);
end