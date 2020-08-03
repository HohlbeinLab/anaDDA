%% Script to make figures 1E of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure1E.mat']); %Loads input file if stored in same location as this script
koffrange = [1 10 100];
konrange = [1 10 100];
for i = koffrange
    for j = konrange
        input.koff1_A = i;
        input.kon1_A = j;
        Comparesimulationwiththeory(input)
    end
end