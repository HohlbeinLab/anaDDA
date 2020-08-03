%% Script to make figures 2A-C of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure2ABC.mat']); %Loads input file if stored in same location as this script
koffrange = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];         % Range of koff values  0.1 0.2 0.5 1 2 5                                        % Range of koff values
Drange = 4;                                                        % Range of Dfree values
Nparticlerange = [1000 5000 10000 50000];                          % Range of Nparticles values
[koffout,konout,Dfreeout,fraction] = Varyparameterssimulationcomparison(input,koffrange,Drange,Nparticlerange);
figure
hold on
plotgeometricvariance(koffout(:,1,:),koffrange,koffrange)
plotgeometricvariance(koffout(:,2,:),koffrange,koffrange)
plotgeometricvariance(koffout(:,3,:),koffrange,koffrange)
plotgeometricvariance(koffout(:,4,:),koffrange,koffrange)
figure
hold on
plotgeometricvariance(konout(:,1,:),koffrange,koffrange)
plotgeometricvariance(konout(:,2,:),koffrange,koffrange)
plotgeometricvariance(konout(:,3,:),koffrange,koffrange)
plotgeometricvariance(konout(:,4,:),koffrange,koffrange)
figure
hold on
plotgeometricvariance(Dfreeout(:,1,:),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,2,:),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,3,:),Drange,koffrange)
plotgeometricvariance(Dfreeout(:,4,:),Drange,koffrange)

function [koffout,konout,Dfreeout,fraction] = Varyparameterssimulationcomparison(input,koffrange,Drange,Nparticlerange)
%% This function is made to test anaDDA over a range of parameters with simulation
% Requires input file that can me made with the script Generateinputfile
% and adjusted to match desired input parameters 

repeats = 5;                                                       % # times repeated to generate average and std estimates

for l = 1:numel(Nparticlerange)
    input.Nparticles = Nparticlerange(l);
    for j = 1:numel(Drange)
        input.Dfree_A = Drange(j);
        maxDfree = input.upperDfree+input.sigmaerror^2/min(input.frametimerange);
        maxrangeD =-log(maxDfree*1e-10)*maxDfree;
        rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;

        for z = 1:numel(input.frametimerange)
            locerror = input.sigmaerror.^2/input.frametimerange(z);
            [input.dist(z).locerrorpdf,input.dist(z).locerrorpdfcorrected] = makelocerrordistributions(rangeD,locerror,input);
        end


        for i = 1:numel(koffrange)
            input.koff1_A = koffrange(i);
            input.kon1_A = koffrange(i);

            for k = 1:repeats
                [parameters,~,~,KSSTAT] = Comparesimulationwiththeory(input)
                koffout(i,j,k,l) = parameters(2);
                konout(i,j,k,l) = parameters(3);
                Dfreeout(i,j,k,l)= parameters(4);
                if input.numberofspecies >1
                fraction(i,j,k,l)= parameters(1);
                else
                fraction(i,j,k,l)= NaN; 
                end
            end
        end
    end
end

end