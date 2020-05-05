function [parameters,bootstrapparamstd,KSSTAT] = anaDDA(input,inputdata)
%%JH_testBranch!!
%% Basic starting function for anaDDA
% anaDDA requires an input of parameters and data. This function is there
% to supply this set of parameters and data in different ways. 
% Either the user can add an input file that they previously generated with 
% Generateinputfile function (also required for more advanced input changes)
% and a dataset of D values or tracks (described in manual or can see the example files).
% Otherwise the user can run this
% code and prompts will ask the user to fill in the most important
% parameters and a link to the dataset they want to analyze.

% Input:
% - input file with all parameters
% - inputdata file with either tracks

% Output:
% - parameters: Kinetic parameters fo fit (fraction, koff, kon, Dfree)
% - bootstrapparamstd: Standard deviation of found parameters with bootstrapping
% - KSSTAT: Kolmogorov Statistic comparing the distribution of data to the distribution of the theoretical predicted anaDDA 
% - Plots of fitted D distributions

%% Parameters prompt (only if users don't supply their own input file) 
if ~(exist('input') == 1)         % Checks if users already supply an input file
    input = Generateinputfile;  %calls function
    opts.Interpreter = 'tex';
    answers = inputdlg({'(1) How many species do you want to fit?',... %JH could be also non-interconverting?!
        '(2) What is the maximum expected diffusion coefficient \it{D} \rm(µm^2/s)?',...
        '(3) What is the estimated localization error in your measurement (µm)?',...
        '(4) Are there any (spherical/rod-shaped) confinement boundaries (y/n)?',...
        '(5) Are there any maximum tracking windows applied (y/n)?',...
        '(6) Do you want to fix any parameters in the fit (y/n)?'},...
        'Input parameters', 1, {'1','10','0.03','n','n','n'}, opts);
    input.numberofspecies = str2num(answers{1});
    input.upperDfree = str2num(answers{2});
    input.sigmaerror = str2num(answers{3});
    
    if answers{4} == 'y'
        input.confinement = 1;
        answersconfinement = inputdlg({'What is the radius of confinement (µm)?',...
            'What is the length of confinement (if cell is spherical then same as above)?'},...
            'Input confinement',1,{'0.5','3'}, opts);
        input.radiusofcell = str2num(answersconfinement{1});
        input.lengthcell = str2num(answersconfinement{2});
    elseif answers{4} == 'n'
        input.confinement = 0;
        input.radiusofcell = 3000; %careful: hard-coded number
        input.lengthcell = 3000; %careful: hard-coded number
    else 
        error('you have to reply with y or n')
    end
    if answers{5} == 'y'
        input.compensatetracking = 1;
        answertrack = inputdlg('What is the length of the tracking window (µm)?',...
            'Input tracking window',1,{'1'});
        input.trackingwindow = str2num(answertrack{1}); 
    else
        input.compensatetracking = 0;
        input.trackingwindow = 3000; %careful: hard-coded number
    end
    if answers{6} == 'y'
        for i = 1:input.numberofspecies
            answerfixspecies1 = inputdlg(['Confinment parameters species ' num2str(i) ...
                ' (fraction koff (s^-1) kon (s^-1) Dfree (µm^2/s)); (-1 if not fixed)?'] ,...
                ['Input fixed parameters species ' num2str(i)],1,{'-1 20 20 4'});
            input.fixedparameters(i,:) = str2num(answerfixspecies1{1});
            if input.fixedparameters(i,4) == 0
                input.fixedparameters(i,2) = 0.00001;
                input.fixedparameters(i,3) = 100000;
                input.fixedparameters(i,4) = 1;
            end
            
        end
    else
        input.fixedparameters(:) = -1;
    end
end

%% Prompt to select input data
if ~exist('inputdata')
[tracksfilename, tracksPathname] = uigetfile('MatFile data:','MultiSelect', 'on');
inputdata = importdata([tracksPathname tracksfilename]);
end


%% Change tracks to list of diffusion coefficients
%Structure of tracks (2D) 
%   column 1: x position (um)
%   column 2: y position (um)
%   column 3: frame number
%   column 4: track id
%   column 5: frame time (s)

%Structure of D 
%   column 1: avg. diff.coef. (um2/s)
%   column 2: tracklength (frames)
%   column 3: frame time (s)

dataformat = size(inputdata);
if dataformat(2) == 5
    D = [];
    input.frametimerange = unique(inputdata(:,5));
    for j = input.frametimerange
    tracks = inputdata(inputdata(:,5)==j,:);
    input.frametime = j;
    truncation = 8; %careful: hard-coded number!
    [Dtemp] = GenerateDfromtracks(tracks,input,truncation);
    D = [D Dtemp];
    end
else
    D = inputdata';
    input.frametimerange = unique(D(3,:));
end

%% Pre-calculation of distributions for localization error (to speed up later fitting)
for j = 1:numel(input.frametimerange)
maxDfree = input.upperDfree+input.sigmaerror^2/min(input.frametimerange(j));
maxrangeD =-log(maxDfree*1e-10)*maxDfree;
rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;

locerror = input.sigmaerror.^2/input.frametimerange(j);
input.frametime = input.frametimerange(j);
[input.dist(j).locerrorpdf,input.dist(j).locerrorpdfcorrected] = makelocerrordistributions(rangeD,locerror,input);
end

%% Fitting of data with MLE and bootstrapping
[parameters,bootstrapparamstd] = MLEanaDDA(D,rangeD,input);

if numel(input.framerange) == 8 && input.plotlog == true
    f = figure;
    f.WindowState = 'maximized';
end
%% Plotting functions and calculation of KSSTAT
for i = 1:numel(input.frametimerange)
    for j = input.framerange
        framenr = j;
        if exist('tracks')
            truncation = framenr;
            [D] = GenerateDfromtracks(tracks,input,truncation);
        end
        if input.plotlog == true
        if numel(input.framerange) == 8
             subplot(4,2,j,'Parent',f)
             title(['D distribution for track length ' num2str(j) ' steps'])
            hold on
            plotlog(framenr,parameters,D(:,D(3,:)==input.frametime), input,rangeD,bootstrapparamstd,1,false)
         else
            figure
            title(['D distribution for track length ' num2str(j) ' steps'])
            hold on
            plotlog(framenr,parameters,D(:,D(3,:)==input.frametime), input,rangeD,bootstrapparamstd,1,true)
        end
        end

        if input.KSstats == true % Only calculate KSstats if true       
            [~,KSSTAT(i,j)]=kstestanaDDA(framenr,parameters,D(:,D(3,:)==input.frametime), input,rangeD,i);
        else
            KSSTAT = 0;
        end
    end
end

% Make output log file
t = datestr(now,'yy_mm_dd_HHMM');
if exist('tracksPathname')
    savefilename = [tracksPathname tracksfilename(1:end-4) '_outputanaDDA_' t '.txt'];
else
    savefilename = ['outputanaDDA_' t '.txt'];
end

fileid = fopen(savefilename,'a');
fprintf(fileid','%s\n\n',['Output file of anaDDA run on ' tracksfilename]);
fprintf(fileid','%s\n', 'Parameters and bootstrap std values for each species');

fprintf(fileid','%s\t%s\t%s\t%s\n',[string('fraction'),string('koff'),string('kon'),string('Dfree')]);
for i = 1:input.numberofspecies
    fprintf(fileid, '%s\n', '----------------------------------------');
    fprintf(fileid','%s\t%s\t%s\t%s\n',string(parameters(i,:)));
    fprintf(fileid','%s\t%s\t%s\t%s\n',string(bootstrapparamstd(i,:)));
    fprintf(fileid, '%s\n', '----------------------------------------');
end
if input.KSstats == true
    fprintf(fileid, '\n%s\n','KSSTATS for the different plots:');
    for i = 1:numel(input.frametimerange)
        for j = input.framerange
            fprintf(fileid, '%s\n',['frametime: ' num2str(input.frametimerange(i)) ' s; track length: ' num2str(input.framerange(j)) ' steps']);
            fprintf(fileid, '%s\n',['KSSTAT: ' num2str(KSSTAT(i,j))]);
        end
    end
end
fclose(fileid);
save([savefilename(1:end-4) '_inputfile'],'input')
profile off
profile viewer
