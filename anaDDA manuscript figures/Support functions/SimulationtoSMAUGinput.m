function [koff,kon,Dfree,D1] = SimulationtoSMAUGinput(input,frametime,k)
totalparticles = input.Nparticles;
distributionstepsize = [28723 20581 14747 10567 7572 5425 3887 8496];
Dlistdata = [];
trajectories={};
trackSMAUG = [];
particlessofar = 0;
for i = input.framerange
    input.NumberofFrames = i;
    input.Nparticles = round(distributionstepsize(i)*totalparticles/100000);
    [~,~,tracks] = SimulationLocalizationandConfinement(input,false);
trackSMAUGtemp = tracks(:,4)+particlessofar;
particlessofar = particlessofar + input.Nparticles;
%trackSMAUGtemp(:,2) = repmat([1:i+1]',input.Nparticles,1);
for j = 1:max(trackSMAUGtemp(:,1))
trackSMAUGtemp((trackSMAUGtemp(:,1)==j),2) = 1:sum(trackSMAUGtemp(:,1)==j);
end  
trackSMAUGtemp(:,3) = 1;
trackSMAUGtemp(:,4) = tracks(:,1);
trackSMAUGtemp(:,5) = tracks(:,2);
trackSMAUG = [trackSMAUG;trackSMAUGtemp];
end
trfile = trackSMAUG;
filepath = mfilename('fullpath');
path = fileparts(filepath);
file =  ['\SimulationSMAUG_Npart_' num2str(totalparticles) '_koff_' num2str(input.koff1_A) '_kon_' num2str(input.kon1_A) '_Dfree_' num2str(input.Dfree_A) '_sigmaerror_' num2str(input.sigmaerror) '_' num2str(k)]; 
save([path, file, '.mat'], 'trfile')
out = SMAUG(trfile,file(2:end));
koff = out.TransMat{end}(2,1).*1/input.frametime;
kon = out.TransMat{end}(1,2).*1/input.frametime;
Dfree = max(out.Dvals{end});
D1 = min(out.Dvals{end});
% res=VB3_HMManalysis2('D:\jnavink\Desktop\Manuscript Analysis\Figure 3\runinput_08_Apr_2019.m',[path, file, '.mat']);
% koff = res.Wbest.est.aMean(1)/frametime;
% kon = res.Wbest.est.aMean(2)/frametime;
% Dfree = res.Wbest.est.DdtMean(2)/frametime - input.sigmaerror^2/frametime;