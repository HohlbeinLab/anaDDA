function [parametersbest,bootstrapparammean,bootstrapparamstd,locerrorparameter,nllbest] = MLEfitDynamic(Dlistdata,input)
%% Used to extract parameters based on MLE estimation. 
% Can estimate parameters for up to three separate species, each with a koff, kon, Dfree (or D1 and D2) and abundance c. 
% Each species adds 4 parameters to be estimated if not restricted. However
% two degrees of freedom are removed due to the following equations:

% Number of degrees of freedom (all the ones that are 0) 
% depend on species and which parameters are fixed
fixedparameters = input.fixedparameters; 
numberofspecies = input.numberofspecies;
fixedparameters(1,1) = 1;
fixedparameterstemp = fixedparameters(1:numberofspecies,:);
indexfittingparameters = fixedparameterstemp==-1;
indexfittingparameters = [indexfittingparameters; zeros(3-numberofspecies,5)];
indexfittingparameters = logical(indexfittingparameters);

% Initialise starting parameters and lower bound and upperbound
%startparameters = [0.3 50 30 4 1;0.3 1 1 4 1; 0.3 1 1 4 1]; 
%JH added
startparameters = input.startparameters;
startparam = startparameters(indexfittingparameters);
lowerbound = zeros(length(startparam),1)+0.000001;
%JH: 07.04.2022
%upperbound = [ones(sum(fixedparameterstemp(:,1)==-1),1);input.upperstartkon*input.upperstartkoff*ones(sum(sum(indexfittingparameters(:,2:3)==1)),1);input.upperDfree*ones(sum(indexfittingparameters(:,4)==1),1);input.upperDfree*ones(sum(indexfittingparameters(:,5)==1),1)];
upperbound = [ones(sum(fixedparameterstemp(:,1)==-1),1);input.upperstartratio*input.upperstartkoff*ones(sum(sum(indexfittingparameters(:,2:3)==1)),1);input.upperDfree*ones(sum(indexfittingparameters(:,4)==1),1);input.upperDfree*ones(sum(indexfittingparameters(:,5)==1),1)];

%% Make lists to keep track what measurement and track length D values come from
Frametimelist = Dlistdata(3,:);
Numberofframes = Dlistdata(2,:);
Dlistdata = Dlistdata(1,:);
[Numberofframes, sortind]=sort(Numberofframes);
Dlistdata = Dlistdata(sortind);
Frametimelist = Frametimelist(sortind);
[Frametimelist, sortind]=sort(Frametimelist);
Dlistdata = Dlistdata(sortind);
Numberofframes = Numberofframes(sortind);

%% Generate frequency table of the lists above and confined function if required
for j = 1:numel(input.frametimerange)
    table = tabulate(Numberofframes(Frametimelist==input.frametimerange(j)));
    frequency(j,:) = table(:,2);
    input.frametime = input.frametimerange(j);
    if input.confinement == true
        [fx(:,j),fy(:,j)] = Generateconfinedfunction(0:0.05:input.upperDfree,input);
    else
        fx = NaN;
        fy = NaN;
    end 
end

    
% Fitting function
custompdf = @(Dlistdata,varargin) pdfDvaluesMLE(Dlistdata, Numberofframes,input, fixedparameters,indexfittingparameters,fx,fy,frequency,varargin);


%% Start values for MLE
mincyclenumber = input.cyclenumber;
lowerstartkoff = input.lowerstartkoff;
upperstartkoff = input.upperstartkoff;
lowerstartlocerror = input.lowerstartlocerror;

lowerstartkon = input.lowerstartratio;
upperstartkon = input.upperstartratio;
upperstartlocerror = input.upperstartlocerror;
maxDfree = input.upperDfree;

i = 1;
pass = 0;
if input.fitlocerror == 1
    lowerbound = [lowerbound; lowerstartlocerror];
    upperbound = [upperbound; upperstartlocerror];
end

while pass == 0
    disp(['Running fitting cycle ' num2str(i) ' of at least ' num2str(mincyclenumber)])
    for j = 1:3
        koffstart(j) = 10^((log10(upperstartkoff)-log10(lowerstartkoff))*rand()+log10(lowerstartkoff));
        %konstart(j) = (upperstartkon-lowerstartkon)*rand+lowerstartkon;
        %konstart(j) = 10^(-1 + 2*rand());
        konstart(j) = 10^((log10(upperstartkon)-log10(lowerstartkon))*rand()+log10(lowerstartkon));
        locerrorstart = lowerstartlocerror + (upperstartlocerror-lowerstartlocerror)*rand();
    end
    startparameters = [rand koffstart(1) konstart(1) maxDfree*rand maxDfree*rand;rand koffstart(2) konstart(2) maxDfree*rand maxDfree*rand; rand koffstart(3) konstart(3) maxDfree*rand maxDfree*rand]; 
    startparam = startparameters(indexfittingparameters);
    if input.fitlocerror == 1
        startparam = [startparam; locerrorstart];
    end
    tableout = maketable(startparam,fixedparameters, indexfittingparameters,numberofspecies,input.fitlocerror);
    disp('Starting parameters are:')
    disp(tableout)
    try
        parameters(:,i) = mle(Dlistdata, 'pdf',custompdf,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound);
        nll(:,i) = -sum(log(custompdf(Dlistdata',parameters(:,i)')));
        disp('Found parameters are:')
        [tableout,totalparam(:,i),locerrorparam] = maketable(parameters(:,i),fixedparameters, indexfittingparameters,numberofspecies,input.fitlocerror);
        disp(tableout)
    catch
         warning('error in fit: no result for this cycle')
         parameters(:,i) = 0;
         nll(:,i) = 1e10;
    end
    if i >= mincyclenumber
        nllrank = sort(nll);
        bestrun = find(nll==nllrank(1));
        secondbestrun = find(nll==nllrank(2));
        parametersbest = parameters(:,bestrun(1));
        parameterssecondbest = parameters(:,secondbestrun(1));
        parametersbest = totalparam(:,bestrun(1));
        parameterssecondbest = totalparam(:,secondbestrun(1));
        if any(abs(parametersbest./parameterssecondbest-1)>0.05)
            pass = 0;
            if i > 20*mincyclenumber
                pass = 1;
                warning('Passed without convergence because reached max cycle number (20 times min cycle number)')
            end
        else
            pass = 1;
        end
    end
    i =  i + 1;    
end

parameters(:,nll==0)=[];
nll(:,nll==0)=[];
nllbest = min(nll);
bestrun = find(nll==nllbest);
parametersbest1 = parameters(:,bestrun(1));
disp('End of run, best parameters were:')
tableout = maketable(parametersbest1,fixedparameters, indexfittingparameters,numberofspecies,input.fitlocerror);
disp(tableout)
%   if nllbest >  -sum(log(custompdf(Dlistdata',input.koff1_A,input.kon1_A,input.Dfree_A)))
%         sprintf('not enough cyclenumbers')   
%   end

%% Add fitted parameters to already fixed parameters
parametersbest = fixedparameters;
startparam = parametersbest1; 
if input.fitlocerror == 1
    locerrorparameter = parametersbest1(end);
    parametersbest1 = parametersbest1(1:end-1);
else
    locerrorparameter = NaN;
end

try
parametersbest(indexfittingparameters) = parametersbest1;
catch
parametersbest(indexfittingparameters) = startparam; 
sprintf('no parameters found')  
end


c = parametersbest(:,1);
koff = parametersbest(:,2);
kon = parametersbest(:,3);
Dfree = parametersbest(:,4);

if input.numberofspecies<3
c(input.numberofspecies+1:end) = 0;
end
c(1) = 1 - c(2) - c(3);
parametersbest(1,1) = c(1);

Dmean2 = c(2)* Dfree(2)*(1-kon(2)/(koff(2)+kon(2)));
Dmean3 = c(3)* Dfree(3)*(1-kon(3)/(koff(3)+kon(3)));


% if fixedparameters(1,4) == 0 % If Dfree is not fixed
% parametersbest(1,4) = (meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*(1-(kon(1)/(koff(1)+kon(1)))));
% elseif fixedparameters(1,3) == 0 % If kon is not fixed
%    a = 1-(meanD - locerror-max(Dmean2,0)-max(Dmean3,0))/(c(1)*Dfree(1));
%     parametersbest(1,3) = a*koff(1)/(1-a);
% else 
% end

parametersbest = parametersbest(1:numberofspecies,:);
parametersbest(:,3) = parametersbest(:,2).*parametersbest(:,3);
paramsize = size(parametersbest);

%% Do bootstrapping
if input.bootstrapping == true
    disp('Initialising Bootstrapping')
    AmountofSubsamples = input.numberofbootstraps;
    if input.fitlocerror == 1
        paramsize(2) = paramsize(2)+1;
    end
    bootstrapparameters=zeros(paramsize(1),paramsize(2),AmountofSubsamples);
        
    for i = 1:AmountofSubsamples
        disp(['Running Bootstrap ' num2str(i) ' of ' num2str(AmountofSubsamples)])
        Dbootstrap = [];
        for j = input.frametimerange
            for k = input.framerange
                index = Frametimelist == j & Numberofframes == k;
                Dlistdatatemp = Dlistdata(index);      
                randomindex = randi(numel(Dlistdatatemp),numel(Dlistdatatemp),1);
                Dbootstrap(index) = Dlistdatatemp(randomindex);
            end
        end
        [bootstraptemp] = mle(Dbootstrap, 'pdf',custompdf,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound);
            
        bootstrapparametersbest = fixedparameters;
        
        if input.fitlocerror ==1
            bootlocerror = bootstraptemp(end);
            bootstrapparametersbest(indexfittingparameters) = bootstraptemp(end-1);
            bootstrapparametersbest(:,6) = bootlocerror;
        else
            bootstrapparametersbest(indexfittingparameters) = bootstraptemp;
        end
        c = bootstrapparametersbest(:,1);
        if input.numberofspecies<3
            c(input.numberofspecies+1:end) = 0;
        end
        c(1) = 1 - c(2) - c(3);
        bootstrapparametersbest(1,1) = c(1);
        bootstrapparametersbest = bootstrapparametersbest(1:numberofspecies,:);
        bootstrapparameters(:,:,i)=bootstrapparametersbest;
    end
    bootstrapparamstd = std(bootstrapparameters,0,3);
    bootstrapparammean = mean(bootstrapparameters,3);
else
    paramsize = size(parametersbest);
    bootstrapparamstd = zeros(paramsize(1),paramsize(2));
    bootstrapparammean= parametersbest;
end

function [output] = pdfDvaluesMLE(x, Numberofframes,input, fixedparameters,indexfittingparameters,fx,fy,frequency,varargin)
parameters = fixedparameters;

if input.fitlocerror == 0
    parameters(indexfittingparameters) = cell2mat(varargin{1});
else 
   param =  cell2mat(varargin{1});
   parameters(indexfittingparameters) = param(1:end-1);
end

%% One degree of freedom less because meanD is linked to fonDNA and Dfree of all species
c = parameters(:,1);
koff = parameters(:,2);
kon = parameters(:,3).*koff;
Dfree = max(parameters(:,4),parameters(:,5));
D1 = min(parameters(:,4),parameters(:,5));
if input.numberofspecies<3
c(input.numberofspecies+1:end) = 0;
end
c(1) = 1 - c(2) - c(3);

maxdata = max(x'.*Numberofframes)+1;
maxrangeD = max(maxdata,-log(1e-22)*max(Dfree));
rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;
maxindex = numel(rangeD);
output = zeros(numel(Numberofframes),1);
dataind = 1;

%% Generate tracking window limit 
if input.compensatetracking == true
      maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
      maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end

%% Generation of distribution for each frame time
for j = 1:numel(input.frametimerange)
    input.frametime = input.frametimerange(j);
    
    if input.fitlocerror == 0
        locerror = input.sigmaerror.^2/input.frametime;
    else
        locerror = param(end).^2/input.frametime;
    end
    
    combinedpdf = zeros(maxindex,8); 
    for i = 1:input.numberofspecies
        [pdf]= DDistributiongenerator(koff(i),kon(i),Dfree(i),D1(i),rangeD,locerror,fx,fy,maxDindtracking,input,j);
        pdf = pdf./sum(pdf);
        combinedpdf = combinedpdf + c(i)*pdf;
    end
%% Interpolation of x data with found distribution for each track length
for i = input.framerange
    numberofdata = frequency(j,i);
    ind = i*(x(dataind:dataind+numberofdata-1)-rangeD(1))/(rangeD(2)-rangeD(1))+1;
    outputtemp = interp1(combinedpdf(1:round(max(ind))+100,i),ind,'spline')./(rangeD(2)-rangeD(1));    
    output(dataind:dataind+numberofdata-1) = outputtemp;
    dataind = dataind+numberofdata;
end
end
output = max(1e-99,output);
