%% This function allows for generation of a new locdisttable.m file, which is located in the main anaDDA folder. 
% This table is used to generate the distribution of immobile particles
% that remain bound. Because the 'observed' motion of these particles is
% solely based on the localization error, they are highly correlated and
% therefore need to be adjusted for. This function generates these
% distributions. 

function [rangex,locdist] = generatePhi2
corrskew = 0; % This variable adds additional correlation (potentially brightness related) on top of 1/4 standard correlation of localization error
for n = 2:8
    binsize = 0.0025*n;
    rangex(:,n) = 0:binsize:4000*binsize;
    C = [1 sqrt(1/4+3/4*corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew);sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew);sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew); sqrt(corrskew) sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew);sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew) sqrt(corrskew) sqrt(corrskew);sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew) sqrt(corrskew); sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1 sqrt(1/4+3/4*corrskew); sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(corrskew) sqrt(1/4+3/4*corrskew) 1]; % The correlation matrix
    C = C(1:n,1:n); % The correlation matrix only for numbers of frames 
    eigenvalues = eig(C);
    xvector = -rangex(:,n)'./eigenvalues;
    if n ==2
        N = round(200000/n);
    else
        N = round(100000/n); 
    end
    locdist(:,n) = rangex(:,n)'.^(n-1)/((n+1)/(2^n)*gamma(n)).*Phi2(ones(n,1),n,xvector',N)'*(rangex(2,n)-rangex(1,n));
    %trial(:,n) = trial(:,n)./sum(trial(:,n));
end