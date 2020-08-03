%% Script to make figures S4 of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigureS4.mat']); %Loads input file if stored in same location as this script
%% Variable rod
totalparticles = input.Nparticles;
Dlistdata = [];
celllengths = [3 4 5];
totalparticlescell = [totalparticles*0.2 totalparticles*0.6 totalparticles*0.2];
input.radiusofcell = sqrt(5*4*0.01);
for v = 1:numel(celllengths)
    input.lengthcell = celllengths(v);
    for w = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(w);
        for u = input.framerange
            input.NumberofFrames = u;
            input.Nparticles = round(input.distributionNparticles(u)*totalparticlescell(v));
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = u;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
    end
end
input.lengthcell = 4;
parameters = anaDDA(input,Dlistdata');
%% Uniform rod
Dlistdata = [];
celllengths = [4 4 4];
totalparticlescell = [totalparticles*0.2 totalparticles*0.6 totalparticles*0.2];
input.radiusofcell = sqrt(5*4*0.01);
for v = 1:numel(celllengths)
    input.lengthcell = celllengths(v);
    for w = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(w);
        for u = input.framerange
            input.NumberofFrames = u;
            input.Nparticles = round(input.distributionNparticles(u)*totalparticlescell(v));
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = u;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
    end
end
input.lengthcell = 4;
parameters = anaDDA(input,Dlistdata');
%% Variable sphere
Dlistdata = [];
input.radiusofcell = sqrt(5*4*0.01);
celllengths = [2*0.75*input.radiusofcell 2*input.radiusofcell 2*1.25*input.radiusofcell];
totalparticlescell = [totalparticles*0.2 totalparticles*0.6 totalparticles*0.2];
for v = 1:numel(celllengths)
    input.lengthcell = celllengths(v);
    input.radiusofcell = input.lengthcell/2;
    for w = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(w);
        for u = input.framerange
            input.NumberofFrames = u;
            input.Nparticles = round(input.distributionNparticles(u)*totalparticlescell(v));
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = u;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
    end
end
input.radiusofcell = sqrt(5*4*0.01);
input.lengthcell = 2*input.radiusofcell;
parameters = anaDDA(input,Dlistdata');
%% Uniform sphere
Dlistdata = [];
input.radiusofcell = sqrt(5*4*0.01);
celllengths = [2*input.radiusofcell 2*input.radiusofcell 2*input.radiusofcell];
totalparticlescell = [totalparticles*0.2 totalparticles*0.6 totalparticles*0.2];
input.radiusofcell = sqrt(5*4*0.01);
for v = 1:numel(celllengths)
    input.lengthcell = celllengths(v);
    input.radiusofcell = input.lengthcell/2;
    for w = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(w);
        for u = input.framerange
            input.NumberofFrames = u;
            input.Nparticles = round(input.distributionNparticles(u)*totalparticlescell(v));
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = u;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
    end
end
input.radiusofcell = sqrt(5*4*0.01);
input.lengthcell = 2*input.radiusofcell;
parameters = anaDDA(input,Dlistdata');
