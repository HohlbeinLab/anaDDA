function [koff,kon,Dfree,Numberofstates,D1] = SimulationtoHMMinput(input,frametime,k)
totalparticles = input.Nparticles;
distributionstepsize = [28723 20581 14747 10567 7572 5425 3887 8496];
Dlistdata = [];
trajectories={};
for i = input.framerange
    input.NumberofFrames = i;
    input.Nparticles = round(distributionstepsize(i)*totalparticles/100000);
    [~,~,tracks] = SimulationLocalizationandConfinement(input,false);

    for i=1:max(tracks(:,4))
        rows = (input.NumberofFrames+1)*(i-1)+1:(input.NumberofFrames+1)*(i-1)+1+input.NumberofFrames;
        if numel(rows)>=input.NumberofFrames
            %trajectories{end+1}=tracks(rows,[1 2]);
            trajectories{end+1}=tracks(tracks(:,4)==i,[1 2]);
        end
    end
end
path = 'D:\jnavink\Desktop\Manuscript Analysis\Figure 3';
file =  ['\SimulationvbSPT_Npart_' num2str(totalparticles) '_koff_' num2str(input.koff1_A) '_kon_' num2str(input.kon1_A) '_Dfree_' num2str(input.Dfree_A) '_sigmaerror_' num2str(input.sigmaerror) '_' num2str(k)]; 
save([path, file, '.mat'], 'trajectories')
res=VB3_HMManalysis2('D:\jnavink\Desktop\Manuscript Analysis\Figure 3\runinput_08_Apr_2019.m',[path, file, '.mat']);
%res=VB3_HMManalysis2('D:\jnavink\Desktop\Manuscript Analysis\Figure 3\runinput_04_Oct_2019_5statesmin8.m',[path, file, '.mat']);
Numberofstates = res.Wbest.N;
koff = res.Wbest.est.aMean(1)/frametime;
kon = res.Wbest.est.aMean(2)/frametime;
Dfree = res.Wbest.est.DdtMean(2)/frametime - input.sigmaerror^2/frametime;
D1 = res.Wbest.est.DdtMean(1)/frametime - input.sigmaerror^2/frametime;