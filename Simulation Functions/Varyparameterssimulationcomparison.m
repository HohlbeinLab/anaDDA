function [koffout,konout,Dfreeout,D1out,fraction] = Varyparameterssimulationcomparison(input)
%% This function is made to test anaDDA over a range of parameters with simulation
% Requires input file that can me made with the script Generateinputfile
% and adjusted the range of kinetic parameters below to match desired input parameters 
koffrange = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];                                            % Range of koff values
konrange = 1;                                                   % Range of kon values (kon*koff)
Drange = 4;                                                     % Range of Dfree values
Nparticlerange = [50000];                                       % Range of Nparticles values
repeats = 5;                                                    % # times repeated to generate average and std estimates


for l = 1:numel(Nparticlerange)
    input.Nparticles = Nparticlerange(l);
    for k = 1:numel(Drange)
        input.Dfree_A = Drange(k);
        maxDfree = input.upperDfree+input.sigmaerror^2/min(input.frametimerange);
        maxrangeD =-log(maxDfree*1e-10)*maxDfree;
        rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;


        for i = 1:numel(koffrange)
            input.koff1_A = koffrange(i);
            for j = 1:numel(konrange)
            input.kon1_A = koffrange(i)*konrange(j);

                for m = 1:repeats
                    [parameters,~,~,KSSTAT] = Comparesimulationwiththeory(input)
                    koffout(i,j,k,l,m) = parameters(2);
                    konout(i,j,k,l,m) = parameters(3);
                    Dfreeout(i,j,k,l,m)= parameters(4);
                    D1out(i,j,k,l,m) = parameters(5);
                    if input.numberofspecies >1
                    fraction(i,j,k,l,m)= parameters(1); 
                    else
                    fraction(i,j,k,l,m)= NaN;
                    end
                end
            end
        end
    end
end