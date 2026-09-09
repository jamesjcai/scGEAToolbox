function [X] = sc_simudata(numgenes, numcells, methodtype)
    % SC_SIMUDATA  Simulate single-cell RNA-seq count data
    %
    %   X = sc_simudata(numgenes, numcells, methodtype)
    %
    %   Inputs:
    %     numgenes   - number of genes (default: 500)
    %     numcells   - number of cells (default: 1000)
    %     methodtype - 'simple' (default) or 'lun'
    %
    %   Output:
    %     X - a numgenes x numcells simulated count matrix
    %
    % References:
    %   * Simple NB sampling: Gamma-NB model (Splat/splatter)
    %   * Lun method: pooling across cells (PMID: 27122128)

% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1305-0
if nargin < 1, numgenes = 500; end
if nargin < 2, numcells = 1000; end
if nargin < 3, methodtype = 'simple'; end

% The NB dispersion is set for each gene at phi_i = 0.1. These parameter
% values were chosen to recapitulate aspects of real data; see
% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7
phi = 0.1;


%simData <- function(counts, models = c("Splat", "SplatDrop", "Simple", "Lun",
%                                       "Lun2", "Lun2ZINB", "scDD", "BASiCS"),

switch lower(methodtype)
    case 'simple'
        % Gene means from a gamma, then NB counts about those means.
        % Mean shape  a  Shape parameter for the mean gene expression gamma
        % Mean rate   b  Rate parameter for the mean gene expression gamma
        shape = 2;
        scale = 0.5;
        genemean = gamrnd(shape, scale, [numgenes, 1]);
        X = i_nbcounts(repmat(genemean, 1, numcells), phi);
    case 'lun'
        % PMID:27122128
        % Pooling across cells to normalize single-cell RNA sequencing data
        % with many zero counts
        % https://github.com/Oshlack/splatter/blob/master/R/lun-simulate.R
        shape = 2;
        scale = 0.5;
        genemean = gamrnd(shape, scale, [numgenes, 1]);

        % cell.facs <- 2 ^ rnorm(nCells, sd = 0.5), i.e. the standard
        % deviation belongs INSIDE the exponent. This was
        % 2.^randn(1, numcells)./2, which divides after exponentiating: the
        % factors then had median 0.501 instead of 1.001 and an sd on the
        % log2 scale of 1.001 instead of 0.500, so every simulated library
        % was about half the intended size and the spread of library sizes
        % was twice as wide as the method specifies.
        cellfactor = 2.^(randn(1, numcells)*0.5);

        X = i_nbcounts(genemean*cellfactor, phi);
    otherwise
        error('sc_simudata:unknownMethod', ...
            'Unknown methodtype ''%s''. Use ''simple'' or ''lun''.', ...
            methodtype);
end
end

function X = i_nbcounts(mu, phi)
%I_NBCOUNTS  NB counts with mean MU and dispersion PHI.
%
%   Splatter draws counts as rnbinom(mu = mu, size = 1/phi), which has
%   mean mu and variance mu + phi*mu^2. MATLAB's NBINRND(R, P) is
%   parameterised by a number of successes R and a success probability P,
%   with mean R*(1 - P)/P, so the translation is R = 1/phi and
%   P = R/(R + mu).
%
%   This used to be called as nbinrnd(mu, phi): the gene mean was passed
%   as the SIZE parameter and the dispersion as the SUCCESS PROBABILITY.
%   Both roles were wrong, and the result was not the requested
%   distribution in either its mean or its variance. Measured over 400
%   genes and 3000 cells with gene means of 1: the 'simple' branch
%   returned means of 9.41 -- the mean of nbinrnd(mu, 0.1) is 9*mu -- with
%   an implied dispersion of 1.90 rather than 0.10, and a variance/mean
%   ratio pinned at 1/p = 10 for every gene, which is a fixed
%   overdispersed-Poisson shape rather than an NB with the stated phi. The
%   'lun' branch returned means of 5.74 and an implied dispersion of 3.58.

r = 1/phi;
X = nbinrnd(r, r./(r + mu));
end
