function [tableout,totalparam,locerrorparam] = maketable(fittingparameters,fixedparameters,indexfittingparameters,numberofspecies,locerrorfit)
if locerrorfit == 1
    locerrorparam = fittingparameters(end);
    fittingparameters = fittingparameters(1:end-1);
else
    locerrorparam = NaN;
end
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
else
    if locerrorfit == 1
        locerror = locerrorparam*ones(length(species),1);
        tableout = table(species,fraction,koff,kon,Dfree,locerror);
    else
        tableout = table(species,fraction,koff,kon,Dfree);
    end
    
end
totalparam = totalparam(:);