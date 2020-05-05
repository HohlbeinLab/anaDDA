function [parameters,bootstrapparamstd,Dlistdata,KSSTAT] = Comparesimulationwiththeory(rangeD,Dfit,simulation,input)
% Decides whether you want to use Simulated values (1) or data (~1), which
% should be supplied as input variable
totalparticles = input.Nparticles;
Dlistdata = [];

%% Simulation part
if simulation == true
    for j = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(j);                                               % Range of different frametimes possible
        for i = input.framerange
            input.NumberofFrames = i;
            input.Nparticles = round(input.distributionNparticles(i)*totalparticles);
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = i;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
    end
else
    Dlistdata = Dfit;
end


%% Allow tracking window 
if input.compensatetracking == false
input.trackingwindow = 300;                                                                     % Setting tracking window above 100 cancels taking this into effect                                                     
end

%% Fitting part
if input.nofit == 0
[parameters,bootstrapparamstd] = MLEanaDDA(Dlistdata,rangeD,input);                             % Actual fitting of data
else
parameters = [1 input.koff1_A input.kon1_A input.Dfree_A];                                      % Using same input parameters for both simulation and theoretical distribution
bootstrapparamstd = [0 0 0 0];
end

%% Generate plots and calculate KSSTAT values
for i = 1:numel(input.frametimerange)
    input.frametime = input.frametimerange(i);  
for j = input.framerange
    framenr = j;
    if input.KSstats == true % Only calculate KSstats if true       
        [~,KSSTAT(i,j)]=kstestanaDDA(framenr,parameters,Dlistdata(:,Dlistdata(3,:)==input.frametime), input,rangeD,i);
    else
        KSSTAT = 0;
    end
    if input.plotlog == true % Only plot if this is true
        figure
        hold on
        plotlog(framenr,parameters,Dlistdata(:,Dlistdata(3,:)==input.frametime), input,rangeD,bootstrapparamstd,i)
   end
end
end
KSSTAT = KSSTAT(KSSTAT>0);