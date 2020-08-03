%% Script to make figures 3BCF of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure3BCF.mat']); %Loads input file if stored in same location as this script
%%Sphere and rods
totalparticles = input.Nparticles;
Dlistdata = [];
ratio = [2 8];
size = [2 5 10];
for i = ratio
for j = size
input.radiusofcell = sqrt(j*4*0.01);
input.lengthcell = i*input.radiusofcell;
input.confinement = 1;
Dlistdata = [];
for w = 1:length(input.frametimerange)
    input.frametime= input.frametimerange(w);
    for u = input.framerange
        input.NumberofFrames = u;
        [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
        Dlistdatatemp(2,:) = u;
        Dlistdatatemp(3,:) = input.frametime;
        Dlistdata = [Dlistdata Dlistdatatemp];
    end
end
parameters = anaDDA(input,Dlistdata');
input.confinement = 0;
parameters = anaDDA(input,Dlistdata');
end
end
%%Trackingwindow
trackingwindow = [5 10 20];
for i = trackingwindow
input.trackingwindow = sqrt(i*4*0.01);
Dlistdata = [];
input.compensatetracking = 1;
for w = 1:length(input.frametimerange)
    input.frametime= input.frametimerange(w);
    for u = input.framerange
        input.NumberofFrames = u;
        [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
        Dlistdatatemp(2,:) = u;
        Dlistdatatemp(3,:) = input.frametime;
        Dlistdata = [Dlistdata Dlistdatatemp];
    end
end
parameters = anaDDA(input,Dlistdata');
input.compensatetracking = 0;
parameters = anaDDA(input,Dlistdata');
end
