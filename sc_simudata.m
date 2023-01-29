function [X]=sc_simudata(numgenes,numcells,methodtype)

% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0
if nargin<1, numgenes = 500; end
if nargin<2, numcells = 1000; end
if nargin<3, methodtype='simple'; end
X=zeros(numgenes, numcells);


%simData <- function(counts, models = c("Splat", "SplatDrop", "Simple", "Lun",
%                                       "Lun2", "Lun2ZINB", "scDD", "BASiCS"),

switch lower(methodtype)
    case 'simple'
        % Define the shape and scale parameters of the gamma distribution
        % Mean shape	α Shape parameter for the mean gene expression gamma distribution
        % Mean rate	β Rate parameter for the mean gene expression gamma distribution
        shape = 2;
        scale = 0.5;
        % Generate random data from the gamma distribution
        rv = gamrnd(shape, scale, [numgenes, 1]);
        for k=1:numgenes
            r = rv(k);
            p = 0.1;
            X(k,:) = nbinrnd(r,p,1,numcells);   % negative binomial distribution
        end
    case 'lun'
        % https://github.com/Oshlack/splatter/blob/master/R/lun-simulate.R
        shape = 2;
        scale = 0.5;
        rv = gamrnd(shape, scale, [numgenes, 1]);
        %rv = 2.^rv;
        %f = 10+randn(1,numcells)./2;
        %cell.facs <- 2 ^ rnorm(nCells, sd = 0.5)
        f=2.^randn(1,numcells)./2;
        F=rv*f;
        p = 0.1;   % The NB dispersion is also set for each gene at φ i =0.1. 
        % These parameter values were chosen to recapitulate aspects of real data [5] 
        % https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7#Abs1
        
        for k=1:numgenes
            for l=1:numcells
                X(k,l) = nbinrnd(F(k,l),p);   % negative binomial distribution
            end
        end
    otherwise
end

