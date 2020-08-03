%% Script to make figures 2D-F of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure2DEF.mat']); %Loads input file if stored in same location as this script
koffrange = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];
Drange = [4];
Nparticlerange = [50000];
for l = 1:numel(Nparticlerange)
input.Nparticles = Nparticlerange(l);
for j = 1:numel(Drange)
input.Dfree_A = Drange(j);
for i = 1:numel(koffrange)
input.koff1_A = koffrange(i);
input.kon1_A = koffrange(i);
for k = 1:5
[koffvbspt(i,j,k,l) konvbspt(i,j,k,l) Dfreevbspt(i,j,k,l),Numberofstates(i,j,k,l)] = SimulationtoHMMinput(input,0.01,k)
[koffSMAUG(i,j,k,l) konSMAUG(i,j,k,l) DfreeSMAUG(i,j,k,l)] = SimulationtoSMAUGinput(input,0.01,k)
end
end
end
end

figure
hold on
plotgeometricvariance(koffvbspt(:,1,:),koffrange,koffrange)
plotgeometricvariance(koffSMAUG(:,1,:),koffrange,koffrange)
figure
hold on
plotgeometricvariance(konvbspt(:,1,:),koffrange,koffrange)
plotgeometricvariance(konSMAUG(:,1,:),koffrange,koffrange)
figure
hold on
plotgeometricvariance(Dfreevbspt(:,1,:),Drange,koffrange)
plotgeometricvariance(DfreeSMAUG(:,1,:),Drange,koffrange)