%% Script to make figures S1 of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigureS1.mat']); %Loads input file if stored in same location as this script
Comparesimulationwiththeory(input)
input.framerange = 4;
sigmarange = [0.02 0.05];
frametimerange = [0.005 0.01 0.02 0.05];
for i = sigmarange
    for j = frametimerange
        input.sigmaerror = i;
        input.frametimerange = j;
        Comparesimulationwiththeory(input)
    end
end