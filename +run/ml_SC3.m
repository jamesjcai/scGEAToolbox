function [c, cls] = ml_SC3(X, k, plotit, donorm, dolog1)
% SC3 - consensus clustering of single-cell RNA-seq data
% Ref: https://www.nature.com/articles/nmeth.4236

if nargin < 5, dolog1 = true; end
if nargin < 4, donorm = true; end
if nargin < 3
    plotit = false;
end
if nargin < 2
    % [optimk]=fun_num_cluster(X,'type','simlr');
    % [optimk]=fun_num_cluster(X,'type','sc3');
    [optimk] = pkg.e_numclusters(X);
else
    optimk = k;
end

if donorm, [X] = sc_norm(X); end
if dolog1, [X] = log1p(X); end


drange = get_drange(X);
disp('Processing Euclidean distance matrix... 1/3')
Dis = squareform(pdist(X'));
[cls1] = get_clusterarray(Dis, optimk, drange);

disp('Processing Spearman distance matrix... 2/3')
Dis = 1 - corr(X, 'type', 's');
[cls2] = get_clusterarray(Dis, optimk, drange);

disp('Processing Pearson distance matrix... 3/3')
Dis = 1 - corr(X, 'type', 'p');
[cls3] = get_clusterarray(Dis, optimk, drange);

oldpath = pwd;
cleanupCwd = onCleanup(@() cd(oldpath));
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(fileparts(pw1), 'external', 'ml_SC3', 'ClusterPack');
if ~(ismcc || isdeployed), addpath(pth); end

% ClusterPack shells out. SGRAPH writes a partgraph<id>.bat beside the
% graph file and runs it with system(), and that .bat in turn calls
% pmetis or shmetis by bare name. Both rely on the shell resolving a
% bare name from the current directory, which Windows does not do when
% NoDefaultCurrentDirectoryInExePath is set. When it cannot, every
% partition fails with "sgraph: partitioning not successful due to
% external error", each of CSPA, HGPA and MCLA returns a single cluster
% covering every cell, and CLUSTERENSEMBLE picks the best of the three
% at NMI 0 -- which ML_SC3 then returned as its answer, with no error
% and no output argument to signal it, on data whose 24 base
% clusterings had each recovered the truth exactly.
%
% Putting the folder on the process PATH makes the .bat and the two
% executables resolvable by name. Measured on a 120-cell three-group
% fixture: before, one cluster of 120; after, three of 40/40/40, with
% cspa and mcla at NMI 0.983.
oldEnvPath = getenv('PATH');
cleanupEnv = onCleanup(@() setenv('PATH', oldEnvPath));
setenv('PATH', [pth, pathsep, oldEnvPath]);

cd(pth);
cls = [cls1; cls2; cls3];
c = clusterensemble(cls, optimk);
c = c(:);

% Backstop for the case the partitioner still cannot run. One cluster,
% when every base clustering was asked for OPTIMK and found it, is not
% an answer -- and it is the shape a caller is least likely to check.
if optimk > 1 && numel(unique(c)) < 2
    error('run:ml_SC3:consensusFailed', ...
        ['The cluster ensemble collapsed to a single cluster while its ', ...
        '%d base clusterings were each asked for %d. ClusterPack could ', ...
        'not run its graph partitioner; look for "sgraph: partitioning ', ...
        'not successful" above, and check that pmetis and shmetis in ', ...
        '%s are executable.'], size(cls, 1), optimk, pth);
end

if plotit
    % CLUSION documents its arguments as a square SIMILARITY matrix and
    % a cluster label ROW vector, and it checks the second: "if
    % size(cl,1) ~= 1, disp('clusion-error: clustering must be row
    % vector'); return". C was forced to a column one line above, so the
    % check always failed -- plotit printed that third-party line and
    % drew nothing, while ML_SC3 returned normally, so a caller looking
    % only for an error saw success. Measured: 0 figures before and
    % after; passing a row gives 1 figure and 6 handles.
    %
    % The first argument was wrong too. DIS at this point is whichever
    % distance matrix was computed last (Pearson, line 32): 0 on the
    % diagonal, larger where cells are LESS alike, and above 1 for
    % anticorrelated pairs. Feeding that to a similarity plot inverts
    % its block structure. The consensus matrix below is SC3's own
    % object and the thing CLUSION exists to display.
    clusion(i_consensusmatrix(cls), reshape(c, 1, []));
end

end


function S = i_consensusmatrix(cls)
%I_CONSENSUSMATRIX Fraction of the ensemble in which each pair co-clusters.
%   A similarity in [0, 1] with 1 on the diagonal, over the same
%   clusterings CLUSTERENSEMBLE was given.
n = size(cls, 2);
S = zeros(n);
for k = 1:size(cls, 1)
    L = cls(k, :);
    S = S + double(L == L');
end
S = S ./ size(cls, 1);
end


function [cls] = get_clusterarray(Dis, optimk, drange)
[Vs1] = pca(Dis);
[Vs2] = transform_Laplacian(Dis, max(drange));

nD = length(drange);
cls = zeros(2*nD, size(Dis, 2));

textprogressbar('Calculating cluster array: ');
% TEXTPROGRESSBAR keeps its carriage-return state in a persistent, and
% the terminating call below is only reached on success. An interrupted
% or failed run left that state set, and the next ML_SC3 call in the
% session then took the TERMINATION branch on its initialising string --
% after which the first numeric call errored with "The text progress
% must be initialized with a string", in a different function, on a
% later call, with nothing pointing at the run that actually broke.
cleanupBar = onCleanup(@() textprogressbar('', true));

for j = 1:nD
    textprogressbar(100*(j ./ nD));
    idx = kmeans(Vs1(:, 1:drange(j)), optimk, 'MaxIter', 1e9, 'emptyaction', 'singleton', 'replicate', 5);
    cls(2*j-1, :) = idx';
    idx = kmeans(Vs2(:, 1:drange(j)), optimk, 'MaxIter', 1e9, 'emptyaction', 'singleton', 'replicate', 5);
    cls(2*j, :) = idx';
end
textprogressbar('done');
end

function drange = get_drange(X)
% https://www.nature.com/articles/nmeth.4236
% the best clusterings were achieved when d was between 4% and 7% of the number of cells, N (Fig. 1c, Supplementary Fig. 3a and Online Methods).
n = size(X, 2);
drge = round(n.*[0.04, 0.07]);
% 4% of n rounds to 0 below 13 cells, and Vs1(:, 1:0) is empty, so
% kmeans failed with "Expected X to be nonempty" -- a cryptic way to say
% "too few cells". PCA of an n-by-n matrix yields at most n-1
% components and EIGS cannot be asked for n either, so cap there too.
drge = max(1, min(drge, n - 1));
drange = drge(1):drge(2);
% Save and restore the caller's random stream. Seeding the d-range subsample is
% fine; leaving the session parked on that seed is not -- it then
% governs every later tsne, umap and clustering call in the session.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));
rng("shuffle");
if length(drange) > 15
    dx = drange(randperm(length(drange)));
    drange = dx(1:15);
end
end

function [V] = transform_Laplacian(Dis, k)

A = exp(-Dis./max(Dis(:))); % adjacency matrix
%     xD=diag(sum(A).^-0.5);  % D=diag(sum(A)); % d(i) the degree of node i
%     xA=xD*A*xD;             % normalized adjacenty matrix
%     L=eye(size(A,1))-xA;    % also L=xD*(D-A)*xD
[~, L] = i_sbe_laplacian_matrix(A);

% see https://people.orie.cornell.edu/dpw/orie6334/lecture7.pdf
% see https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian_2

%     [V,D]=eig(L);
%     [~,ind]=sort(diag(D));
%     V = V(:,ind);
[V, ~] = eigs(L, k, 'smallestreal');
end


function [L, Lnorm] = i_sbe_laplacian_matrix(A)
% Get graph Laplacian matrix
%
%   L = laplacian(g)
%
% graph Laplacian matrix is defined by L = D - A, where D is vertex degree
% diagonal matrix and A is adjacency matrix.
%
% See also: adjacency
% L = diag(sum(A)) - A;

% Systems Biology & Evolution Toolbox
% Author: James Cai
% Email: jcai@tamu.edu
% Website: https://github.com/jamesjcai/SBEToolbox_lite

% https://github.com/dtuia/KEMA/blob/7378c0fce50a818c2fb59f5de9344ea5c1929fa4/general_routine/laplacian.m
% https://github.com/KavehFathian/clear/blob/e72318ad52082485442f55ad5a7ee9b989b85677/Algorithms/Helpers/NormalizeLap.m
% https://github.com/hungrydoggy/Pinocchio/blob/5664503b210005fa4f1fc053e237ba3ecf6a7945/skeletonizer/matlab/toolbox/compute_mesh_laplacian.m
% https://github.com/fljohnston1/otto-group-product/blob/02a6f35f8c144ed52a2097a1965d561172f2e701/SpectralClustering.m

D = sum(A);
L = diag(D) - A;
if nargout > 1
    D(D ~= 0) = sqrt(1./D(D ~= 0));
    D = diag(D);
    % Lnorm=D*L*D;
    Lnorm = eye(size(A, 1)) - D * A * D; % L = I-D^-1/2*W*D^-1/2
end
end


function textprogressbar(c, resetOnly)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
%
% textprogressbar('', true) clears the persistent state without printing,
% for the error path -- see the onCleanup in GET_CLUSTERARRAY.
% INPUTS:   C   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m
% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version
% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR; %   Carriage return pesistent variable

if nargin > 1 && resetOnly
    strCR = [];
    return
end

% Vizualization parameters
strPercentageLength = 10; %   Length of percentage string (must be >5)
strDotsMaximum = 10; %   The total number of dots in a progress bar

%% Main
if isempty(strCR) && ~ischar(c)
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c)
    % Progress bar - initialization
    fprintf('%s', c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c)
    % Progress bar  - termination
    strCR = [];
    fprintf([c, '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c), '%%'];
    percentageOut = [percentageOut, repmat(' ', 1, strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[', repmat('.', 1, nDots), repmat(' ', 1, strDotsMaximum-nDots), ']'];
    strOut = [percentageOut, dotOut];

    % Print it on the screen
    if strCR == -1
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR, strOut]);
    end

    % Update carriage return
    strCR = repmat('\b', 1, length(strOut)-1);

else
    % Any other unexpected input
    error('Unsupported argument type');
end
end
