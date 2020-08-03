%% Script to make figures S2B of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigureS2B.mat']); %Loads input file if stored in same location as this script
koffrange = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];
Drange = [4];
densityrange = [0 0.01 0.1 0.25];
for l = 1:numel(densityrange)
input.density = densityrange(l)*input.fovsize^2;
for j = 1:numel(Drange)
input.Dfree_A = Drange(j);
for i = 1:numel(koffrange)
input.koff1_A = koffrange(i);
input.kon1_A = koffrange(i);
for k = 1:5
[parameters,~,~,KSSTAT] = Comparesimulationwiththeory(input)
koffout(i,j,k,l) = parameters(2);
konout(i,j,k,l) = parameters(3);
Dfreeout(i,j,k,l)= parameters(4);
D1out(i,j,k,l) = parameters(5);
if input.numberofspecies >1
fraction(i,j,k,l)= parameters(1); 
else
fraction(i,j,k,l)= NaN;
end
if l == 4
[koffvbspt(i,j,k,l) konvbspt(i,j,k,l) Dfreevbspt(i,j,k,l),Numberofstates(i,j,k,l)] = SimulationtoHMMinput(input,0.01,k)
[koffSMAUG(i,j,k,l) konSMAUG(i,j,k,l) DfreeSMAUG(i,j,k,l)] = SimulationtoSMAUGinput(input,0.01,k)
end
end
end
end
end

figure
hold on
plotgeometricvariance(koffout(:,1,:,1),koffrange,koffrange)
plotgeometricvariance(koffout(:,1,:,2),koffrange,koffrange)
plotgeometricvariance(koffout(:,1,:,3),koffrange,koffrange)
plotgeometricvariance(koffout(:,1,:,4),koffrange,koffrange)
plotgeometricvariance(koffvbspt(:,1,:,4),koffrange,koffrange)
plotgeometricvariance(koffSMAUG(:,1,:,4),koffrange,koffrange)
figure
hold on
plotgeometricvariance(konout(:,1,:,1),koffrange,koffrange)
plotgeometricvariance(konout(:,1,:,2),koffrange,koffrange)
plotgeometricvariance(konout(:,1,:,3),koffrange,koffrange)
plotgeometricvariance(konout(:,1,:,4),koffrange,koffrange)
plotgeometricvariance(konvbspt(:,1,:,4),koffrange,koffrange)
plotgeometricvariance(konSMAUG(:,1,:,4),koffrange,koffrange)
figure
hold on
plotgeometricvariance(Dfreeout(:,1,:,1),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,1,:,2),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,1,:,3),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,1,:,4),Drange,koffrange)
plotgeometricvariance(Dfreevbspt(:,1,:,4),Drange,koffrange)
plotgeometricvariance(DfreeSMAUG(:,1,:,4),Drange,koffrange)