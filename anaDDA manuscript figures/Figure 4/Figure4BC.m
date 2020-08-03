%% Script to make figures 4BC of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure4BC.mat']); %Loads input file if stored in same location as this script
%% Simulation part
totalparticles = input.Nparticles;
Dlistdata = [];
maxDfree = input.upperDfree+input.sigmaerror^2/min(input.frametimerange);
maxrangeD =-log(maxDfree*1e-10)*maxDfree;
rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;

%% Simulation part
for j = 1:length(input.frametimerange)
    input.frametime= input.frametimerange(j);                                               % Range of different frametimes possible
    for i = input.framerange
        input.NumberofFrames = i;
        input.Nparticles = round(input.distributionNparticles(i)*totalparticles);
        [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
        Dlistdata = [Dlistdata Dlistdatatemp];
    end
end
D = Dlistdata;
%% Change input 
input.fractionB = 0.2;
input.numberofspecies = 2;
input.koff1_B = 1e-10;
input.kon1_B = 100;
input.koff1_A = 20;
input.kon1_A = 20;
input.koff2_A = 1e10;
input.kon2_A = 1e-10;

[parameters,bootstrapparamstd,KSSTAT] = anaDDA(input,Dlistdata'); 
KSSTAT
input.nofit = false;
input.numberofspecies = 1;
[parameters,bootstrapparamstd,KSSTAT] = anaDDA(input,Dlistdata');
KSSTAT
