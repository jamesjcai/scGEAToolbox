function [c, cls] = mt_SC3(X, k, plotit, donorm, dolog1)
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
    [optimk] = pkg.fun_num_clusters(X);
else
    optimk = k;
end

if donorm, [X] = sc_norm(X); end
if dolog1, [X] = log(X+1); end


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
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, 'external', 'mt_SC3', 'ClusterPack');
if ~(ismcc || isdeployed), addpath(pth); end
cd(pth);
cls = [cls1; cls2; cls3];
c = clusterensemble(cls, optimk);
c = c(:);
cd(oldpath);

if plotit
    clusion(Dis, c);
end

end


function [cls] = get_clusterarray(Dis, optimk, drange)
[Vs1] = pca(Dis);
[Vs2] = transform_Laplacian(Dis, max(drange));
cls = [];
textprogressbar('Calculating cluster array: ');
for j = 1:length(drange)
    textprogressbar(100*(j ./ length(drange)));
    idx = kmeans(Vs1(:, 1:drange(j)), optimk, 'MaxIter', 1e9, 'emptyaction', 'singleton', 'replicate', 5);
    cls = [cls; idx'];
    idx = kmeans(Vs2(:, 1:drange(j)), optimk, 'MaxIter', 1e9, 'emptyaction', 'singleton', 'replicate', 5);
    cls = [cls; idx'];
end
textprogressbar('done');
end

function drange = get_drange(X)
% https://www.nature.com/articles/nmeth.4236
% the best clusterings were achieved when d was between 4% and 7% of the number of cells, N (Fig. 1c, Supplementary Fig. 3a and Online Methods).
n = size(X, 2);
drge = round(n.*[0.04, 0.07]);
drange = drge(1):drge(2);
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


function textprogressbar(c)
% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
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