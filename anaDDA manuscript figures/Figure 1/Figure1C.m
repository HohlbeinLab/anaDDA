%% Script to make figures 1C of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure1C.mat']); %Loads input file if stored in same location as this script
input.Nparticles = 1000;
Comparesimulationwiththeory(input)
input.Nparticles = 10000;
Comparesimulationwiththeory(input)
input.Nparticles = 100000;
Comparesimulationwiththeory(input)