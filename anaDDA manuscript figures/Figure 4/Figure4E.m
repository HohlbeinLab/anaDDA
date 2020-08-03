%% Script to make figures 4E of Vink et al.(2020) Extracting transition rates in particle tracking using analytical diffusion distribution analysis
filepath = mfilename('fullpath');
filepath = fileparts(filepath);
load([filepath '\inputFigure4E.mat']); %Loads input file if stored in same location as this script

frametimeranges = {0.01,0.05,[0.01 0.05]};
nparticle = {100000,100000,50000};

for ii = 1:3
input.frametime = 0.01;
input.frametimerange = frametimeranges{ii};
input.framerange = 4;
input.Nparticles=  nparticle{ii};
Dlistdata = [];

for j = 1:length(input.frametimerange)
        input.frametime= input.frametimerange(j);                                               % Range of different frametimes possible
        for i = input.framerange
            input.NumberofFrames = i;
            [Dlistdatatemp] = SimulationLocalizationandConfinement(input,false);
            Dlistdatatemp(2,:) = i;
            Dlistdatatemp(3,:) = input.frametime;
            Dlistdata = [Dlistdata Dlistdatatemp];
        end
end
koffrange = [10:4:50];
konrange = [10:4:50];
Dfreerange = [0.6:0.08:1.4];
immobilerange = [0:0.02:0.2];
maxDfree = input.upperDfree+input.sigmaerror^2/min(input.frametimerange);
        maxrangeD =-log(maxDfree*1e-10)*maxDfree;
        rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;

        for z = 1:numel(input.frametimerange)
            locerror = input.sigmaerror.^2/input.frametimerange(z);
            [input.dist(z).locerrorpdf,input.dist(z).locerrorpdfcorrected] = makelocerrordistributions(rangeD,locerror,input);
        end
KSSTAT = [];
for i = 1:numel(koffrange)
    i
    for j = 1:numel(konrange)
        for k = 1:numel(Dfreerange)
            for l= 1:numel(immobilerange)
        koff = koffrange(i);
        kon = konrange(j);
        Dfree = Dfreerange(k);    
        immobilefraction = immobilerange(l); 
        parameters = [1-immobilefraction koff kon Dfree 0; immobilefraction 1e-7 100 0.5 0];
        
        [~,KSSTAT(i,j,k,l)]=kstestanaDDA(input.framerange,parameters,Dlistdata, input,rangeD,1);
            end
        end
    end
end
KSSTATresults{ii} = KSSTAT;
KSSTAT
end
for i = 1:3
figure
imagesc(log10(squeeze(min(squeeze(min(KSSTATresults{i},[],2)),[],1))),[-2.5 -1.5])
colormap('gray')
colorbar
ylabel('Dfree')
xlabel('immobilefraction')
xticklabels(strsplit(num2str(immobilerange)))
yticklabels(strsplit(num2str(Dfreerange)))
end
