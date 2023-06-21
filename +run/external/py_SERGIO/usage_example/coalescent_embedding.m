function coords = coalescent_embedding(x, pre_weighting, dim_red, angular_adjustment, dims)

% Authors:
% - main code: Alessandro Muscoloni, 2017-09-21
% - support functions: indicated at the beginning of the function

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, J. M. Thomas, C. V. Cannistraci

% Reference:
% A. Muscoloni, J. M. Thomas, S. Ciucci, G. Bianconi, and C. V. Cannistraci,
% "Machine learning meets complex networks via coalescent embedding in the hyperbolic space",
% Nature Communications 8, 1615 (2017). doi:10.1038/s41467-017-01825-5

% The time complexity of the algorithms is O(N^3).

%%% INPUT %%%
% x - adjacency matrix of the network, which must be:
%   symmetric, zero-diagonal, one connected component, not fully connected;
%   the network can be weighted
%
% pre_weighting - rule for pre-weighting the matrix, the alternatives are:
%   'original' -> the original weights are considered;
%                 NB: they should suggest distances and not similarities
%   'reverse'  -> the original weights reversed are considered;
%                 NB: to use when they suggest similarities
%   'RA1'      -> Repulsion-Attraction v1
%   'RA2'      -> Repulsion-Attraction v2
%   'EBC'      -> Edge-Betweenness-Centrality
%
% dim_red - dimension reduction technique, the alternatives are:
%   'ISO'   -> Isomap (valid for 2D and 3D)
%   'ncISO' -> noncentered Isomap (valid for 2D and 3D)
%   'LE'    -> Laplacian Eigenmaps (valid for 2D and 3D)
%   'MCE'   -> Minimum Curvilinear Embedding (only valid for 2D)
%   'ncMCE' -> noncentered Minimum Curvilinear Embedding (only valid for 2D)
%
% angular_adjustment - method for the angular adjustment, the alternatives are:
%   'original' -> original angular distances are preserved (valid for 2D and 3D)
%   'EA'       -> equidistant adjustment (only valid for 2D)
%   
% dims - dimensions of the hyperbolic embedding space, the alternatives are:
%   2 -> hyperbolic disk
%   3 -> hyperbolic sphere

%%% OUTPUT %%%
% coords - polar or spherical hyperbolic coordinates of the nodes
%   in the hyperbolic disk they are in the form: [theta,r]
%   in the hyperbolic sphere they are in the form: [azimuth,elevation,r]
%   for details see the documentation of the MATLAB functions
%   "cart2pol" and "cart2sph"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
validateattributes(x, {'numeric'}, {'square','finite','nonnegative'});
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
 validateattributes(pre_weighting, {'char'}, {});
validateattributes(dim_red, {'char'}, {});
validateattributes(angular_adjustment, {'char'}, {});
validateattributes(dims, {'numeric'}, {'scalar','integer','>=',2,'<=',3});
if ~any(strcmp(pre_weighting,{'original','reverse','RA1','RA2','EBC'}))
    error('Possible pre-weighting rules: ''original'',''reverse'',''RA1'',''RA2'',''EBC''.');
end
if dims == 2
    if ~any(strcmp(dim_red,{'ISO','ncISO','MCE','ncMCE','LE'}))
        error('Possible dimension reduction techniques in 2D: ''ISO'', ''ncISO'', ''MCE'', ''ncMCE'', ''LE''.');
    end
    if ~any(strcmp(angular_adjustment,{'original','EA'}))
        error('Possible angular adjustment methods in 2D: ''original'', ''EA''.');
    end
elseif dims == 3
    if ~any(strcmp(dim_red,{'ISO','ncISO','LE'}))
        error('Possible dimension reduction techniques in 3D: ''ISO'', ''ncISO'', ''LE''.');
    end
    if ~any(strcmp(angular_adjustment,{'original'}))
        error('Possible angular adjustment methods in 3D: ''original''.');
    end    
end

% pre-weighting
if strcmp(pre_weighting,'original')
    xw = x;
elseif strcmp(pre_weighting,'reverse')
    xw = reverse_weights(x);
elseif strcmp(pre_weighting,'RA1')
    xw = RA1_weighting(double(x>0));
elseif strcmp(pre_weighting,'RA2')
    xw = RA2_weighting(double(x>0));
elseif strcmp(pre_weighting,'EBC')
    xw = EBC_weighting(double(x>0));
end

% dimension reduction and set of hyperbolic coordinates
if dims == 2
    coords = zeros(size(x,1),2);
    if strcmp(dim_red,'ISO')
        coords(:,1) = set_angular_coordinates_ISO_2D(xw, angular_adjustment);
    elseif strcmp(dim_red,'ncISO')
        coords(:,1) = set_angular_coordinates_ncISO_2D(xw, angular_adjustment);
    elseif strcmp(dim_red,'MCE')
        coords(:,1) = set_angular_coordinates_MCE_2D(xw, angular_adjustment);
    elseif strcmp(dim_red,'ncMCE')
        coords(:,1) = set_angular_coordinates_ncMCE_2D(xw, angular_adjustment);
    elseif strcmp(dim_red,'LE')
        coords(:,1) = set_angular_coordinates_LE_2D(xw, angular_adjustment);
    end
    coords(:,2) = set_radial_coordinates(x);
elseif dims == 3
    coords = zeros(size(x,1),3);
    if strcmp(dim_red,'ISO')
        [coords(:,1),coords(:,2)] = set_angular_coordinates_ISO_3D(xw);
    elseif strcmp(dim_red,'ncISO')
        [coords(:,1),coords(:,2)] = set_angular_coordinates_ncISO_3D(xw);
    elseif strcmp(dim_red,'LE')
        [coords(:,1),coords(:,2)] = set_angular_coordinates_LE_3D(xw);
    end
    coords(:,3) = set_radial_coordinates(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Support Functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xrev = reverse_weights(x)

xrev = x;
xrev(xrev>0) = abs(x(x>0) - min(x(x>0)) - max(x(x>0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_RA1 = RA1_weighting(x)

n = size(x,1);
cn = x*x;
deg = full(sum(x,1));
x_RA1 = x .* (repmat(deg,n,1) + repmat(deg',1,n) + ...
    (repmat(deg,n,1) .* repmat(deg',1,n))) ./ (1 + cn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_RA2 = RA2_weighting(x)

n = size(x,1);
cn = x*x;
ext = repmat(sum(x,2),1,n) - cn - 1;
x_RA2 = x .* (1 + ext + ext' + ext.*ext') ./ (1 + cn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_EBC = EBC_weighting(x)

[~,x_EBC] = betweenness_centrality(sparse(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = set_angular_coordinates_ISO_2D(xw, angular_adjustment)

% dimension reduction
dr_coords = isomap_graph_carlo(xw, 2, 'yes');

% from cartesian to polar coordinates
% using dimensions 1 and 2 of embedding
[ang_coords,~] = cart2pol(dr_coords(:,1),dr_coords(:,2));
% change angular range from [-pi,pi] to [0,2pi]
ang_coords = mod(ang_coords + 2*pi, 2*pi);

if strcmp(angular_adjustment,'EA')
    ang_coords = equidistant_adjustment(ang_coords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = set_angular_coordinates_ncISO_2D(xw, angular_adjustment)

% dimension reduction
dr_coords = isomap_graph_carlo(xw, 3, 'no');

% from cartesian to polar coordinates
% using dimensions 2 and 3 of embedding
[ang_coords,~] = cart2pol(dr_coords(:,2),dr_coords(:,3));
% change angular range from [-pi,pi] to [0,2pi]
ang_coords = mod(ang_coords + 2*pi, 2*pi);

if strcmp(angular_adjustment,'EA')
    ang_coords = equidistant_adjustment(ang_coords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = set_angular_coordinates_MCE_2D(xw, angular_adjustment)

% dimension reduction
dr_coords = mce(xw, 1, 'yes');

if strcmp(angular_adjustment,'original')
    % circular adjustment of dimension 1
    ang_coords = circular_adjustment(dr_coords(:,1));
elseif strcmp(angular_adjustment,'EA')
    % equidistant adjustment of dimension 1
    ang_coords = equidistant_adjustment(dr_coords(:,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = set_angular_coordinates_ncMCE_2D(xw, angular_adjustment)

% dimension reduction
dr_coords = mce(xw, 2, 'no');

if strcmp(angular_adjustment,'original')
    % circular adjustment of dimension 2
    ang_coords = circular_adjustment(dr_coords(:,2));
elseif strcmp(angular_adjustment,'EA')
    % equidistant adjustment of dimension 2
    ang_coords = equidistant_adjustment(dr_coords(:,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = set_angular_coordinates_LE_2D(xw, angular_adjustment)

% dimension reduction
st = triu(full(xw),1);
st = mean(st(st>0));
heat_kernel = zeros(size(xw));
heat_kernel(xw>0) = exp(-((xw(xw>0)./st).^2));
dr_coords = leig_graph_carlo_classical(heat_kernel, 2, 'no');

% from cartesian to polar coordinates
% using dimensions 2 and 3 of embedding
% (dimensions 1 and 2 in the code since the first is skipped by the function)
[ang_coords,~] = cart2pol(dr_coords(:,1),dr_coords(:,2));
% change angular range from [-pi,pi] to [0,2pi]
ang_coords = mod(ang_coords + 2*pi, 2*pi);

if strcmp(angular_adjustment,'EA')
    ang_coords = equidistant_adjustment(ang_coords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = equidistant_adjustment(coords)

% sort input coordinates
[~,idx] = sort(coords);
% assign equidistant angular coordinates in [0,2pi[ according to the sorting
angles = linspace(0, 2*pi, length(coords)+1);
ang_coords(idx) = angles(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang_coords = circular_adjustment(coords)

% scale the input coordinates into the range [0,2pi]
n = length(coords);
m = 2*pi*(n-1)/n;
ang_coords = ((coords - min(coords)) ./ (max(coords) - min(coords))) * m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azimuth, elevation] = set_angular_coordinates_ISO_3D(xw)

% dimension reduction
dr_coords = isomap_graph_carlo(xw, 3, 'yes');

% from cartesian to spherical coordinates
% using dimensions 1-3 of embedding
[azimuth,elevation,~] = cart2sph(dr_coords(:,1),dr_coords(:,2),dr_coords(:,3));
% change angular range from [-pi,pi] to [0,2pi]
azimuth = mod(azimuth + 2*pi, 2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azimuth, elevation] = set_angular_coordinates_ncISO_3D(xw)

% dimension reduction
dr_coords = isomap_graph_carlo(xw, 4, 'no');

% from cartesian to spherical coordinates
% using dimensions 2-4 of embedding
[azimuth,elevation,~] = cart2sph(dr_coords(:,2),dr_coords(:,3),dr_coords(:,4));
% change angular range from [-pi,pi] to [0,2pi]
azimuth = mod(azimuth + 2*pi, 2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azimuth, elevation] = set_angular_coordinates_LE_3D(xw)

% dimension reduction
st = triu(full(xw),1);
st = mean(st(st>0));
heat_kernel = zeros(size(xw));
heat_kernel(xw>0) = exp(-((xw(xw>0)./st).^2));
dr_coords = leig_graph_carlo_classical(heat_kernel, 3, 'no');

% from cartesian to spherical coordinates
% using dimensions 2-4 of embedding (the first is skipped by the LE function)
[azimuth,elevation,~] = cart2sph(dr_coords(:,1),dr_coords(:,2),dr_coords(:,3));
% change angular range from [-pi,pi] to [0,2pi]
azimuth = mod(azimuth + 2*pi, 2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radial_coordinates = set_radial_coordinates(x)

n = size(x,1);
deg = full(sum(x>0,1));
if all(deg == deg(1))
   error('All the nodes have the same degree, the degree distribution cannot fit a power-law.'); 
end

% fit power-law degree distribution
gamma_range = 1.01:0.01:10.00;
small_size_limit = 100;
if length(deg) < small_size_limit
    gamma = plfit(deg, 'finite', 'range', gamma_range);
else
    gamma = plfit(deg, 'range', gamma_range);
end
beta = 1 / (gamma - 1);

% sort nodes by decreasing degree
[~,idx] = sort(deg, 'descend');

% for beta > 1 (gamma < 2) some radial coordinates are negative
radial_coordinates = zeros(1, n);
radial_coordinates(idx) = max(0, 2*beta*log(1:n) + 2*(1-beta)*log(n));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [f, time] = leig_graph_carlo_classical(x, d, centring)
% Maps the high-dimensional samples in 'x' to a low dimensional space using
% Laplacian Eigenmaps (coded 27-JANUARY-2013 by Gregorio Alanis-Lobato)

t = tic;

graph = max(x, x');

% Kernel centering
if strcmp(centring, 'yes')
    graph=kernel_centering(graph); %Compute the centred MC-kernel
end

D = sum(graph, 2); %Degree values
D = diag(D); %Degree matrix

% Graph laplacian
L = D - graph;

% Solving the generalised eigenvalue problem L*f = lambda*D*f
[f, ~] = eig(L, D); 
f = real(f(:, 2:d+1));

time = toc(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [s,time] = isomap_graph_carlo(x, n, centring)

%INPUT
%   x => Distance or correlation matrix x
%   n => Dimension into which the data is to be projected
%   centring => 'yes' is x should be centred or 'no' if not
%OUTPUT
%   s => Sample configuration in the space of n dimensions

t = tic;

% initialization
x = max(x, x');

% Iso-kernel computation
% kernel=graphallshortestpaths(sparse(x),'directed','false'); 
kernel=distances(graph(sparse(x))); 

clear x
kernel=max(kernel,kernel'); 

% Kernel centering
if strcmp(centring, 'yes')
    kernel=kernel_centering(kernel); %Compute the centred Iso-kernel
end

% Embedding 
[~,L,V] = svd(kernel, 'econ');
sqrtL = sqrt(L(1:n,1:n)); clear L
V = V(:,1:n);
s = real((sqrtL * V')');

time = toc(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [s, time] = mce(x, n, centring)
%Given a distance or correlation matrix x, it performs Minimum Curvilinear 
%Embedding (MCE) or non-centred MCE (ncMCE) (coded 27-SEPTEMBER-2011 by 
%Carlo Cannistraci)

%INPUT
%   x => Distance or correlation matrix x
%   n => Dimension into which the data is to be projected
%   centring => 'yes' is x should be centred or 'no' if not
%OUTPUT
%   s => Sample configuration in the space of n dimensions

t = tic;

% initialization
x = max(x, x');

% MC-kernel computation
kernel=graphallshortestpaths(graphminspantree(sparse(x),'method','kruskal'),'directed','false'); 

clear x
kernel=max(kernel,kernel'); 

% Kernel centering
if strcmp(centring, 'yes')
    kernel=kernel_centering(kernel); %Compute the centred MC-kernel
end

% Embedding 
[~,L,V] = svd(kernel, 'econ');
sqrtL = sqrt(L(1:n,1:n)); clear L
V = V(:,1:n);
s = (sqrtL * V')';
s=real(s(:,1:n));

time = toc(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = kernel_centering(D)

% 2011-09-27 - Carlo Vittorio Cannistraci

%%% INPUT %%%
% D - Distance matrix

%%% OUTPUT %%%
% D - Centered distance matrix

% Centering
N = size(D,1);
J = eye(N) - (1/N)*ones(N);
D = -0.5*(J*(D.^2)*J);

% Housekeeping
D(isnan(D)) = 0;
D(isinf(D)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, xmin, L]=plfit(x, varargin)
% PLFIT fits a power-law distributional model to data.
%    Source: http://www.santafe.edu/~aaronc/powerlaws/
% 
%    PLFIT(x) estimates x_min and alpha according to the goodness-of-fit
%    based method described in Clauset, Shalizi, Newman (2007). x is a 
%    vector of observations of some quantity to which we wish to fit the 
%    power-law distribution p(x) ~ x^-alpha for x >= xmin.
%    PLFIT automatically detects whether x is composed of real or integer
%    values, and applies the appropriate method. For discrete data, if
%    min(x) > 1000, PLFIT uses the continuous approximation, which is 
%    a reliable in this regime.
%   
%    The fitting procedure works as follows:
%    1) For each possible choice of x_min, we estimate alpha via the 
%       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
%       goodness-of-fit statistic D.
%    2) We then select as our estimate of x_min, the value that gives the
%       minimum value D over all values of x_min.
%
%    Note that this procedure gives no estimate of the uncertainty of the 
%    fitted parameters, nor of the validity of the fit.
%
%    Example:
%       x = (1-rand(10000,1)).^(-1/(2.5-1));
%       [alpha, xmin, L] = plfit(x);
%
%    The output 'alpha' is the maximum likelihood estimate of the scaling
%    exponent, 'xmin' is the estimate of the lower bound of the power-law
%    behavior, and L is the log-likelihood of the data x>=xmin under the
%    fitted power law.
%    
%    For more information, try 'type plfit'
%
%    See also PLVAR, PLPVA

% Version 1.0    (2007 May)
% Version 1.0.2  (2007 September)
% Version 1.0.3  (2007 September)
% Version 1.0.4  (2008 January)
% Version 1.0.5  (2008 March)
% Version 1.0.6  (2008 July)
% Version 1.0.7  (2008 October)
% Version 1.0.8  (2009 February)
% Version 1.0.9  (2009 October)
% Version 1.0.10 (2010 January)
% Version 1.0.11 (2012 January)
% Copyright (C) 2008-2012 Aaron Clauset (Santa Fe Institute)
% Distributed under GPL 2.0
% http://www.gnu.org/copyleft/gpl.html
% PLFIT comes with ABSOLUTELY NO WARRANTY
% 
% Notes:
% 
% 1. In order to implement the integer-based methods in Matlab, the numeric
%    maximization of the log-likelihood function was used. This requires
%    that we specify the range of scaling parameters considered. We set
%    this range to be [1.50 : 0.01 : 3.50] by default. This vector can be
%    set by the user like so,
%    
%       a = plfit(x,'range',[1.001:0.001:5.001]);
%    
% 2. PLFIT can be told to limit the range of values considered as estimates
%    for xmin in three ways. First, it can be instructed to sample these
%    possible values like so,
%    
%       a = plfit(x,'sample',100);
%    
%    which uses 100 uniformly distributed values on the sorted list of
%    unique values in the data set. Second, it can simply omit all
%    candidates above a hard limit, like so
%    
%       a = plfit(x,'limit',3.4);
%    
%    Finally, it can be forced to use a fixed value, like so
%    
%       a = plfit(x,'xmin',3.4);
%    
%    In the case of discrete data, it rounds the limit to the nearest
%    integer.
% 
% 3. When the input sample size is small (e.g., < 100), the continuous 
%    estimator is slightly biased (toward larger values of alpha). To
%    explicitly use an experimental finite-size correction, call PLFIT like
%    so
%    
%       a = plfit(x,'finite');
%    
%    which does a small-size correction to alpha.
%
% 4. For continuous data, PLFIT can return erroneously large estimates of 
%    alpha when xmin is so large that the number of obs x >= xmin is very 
%    small. To prevent this, we can truncate the search over xmin values 
%    before the finite-size bias becomes significant by calling PLFIT as
%    
%       a = plfit(x,'nosmall');
%    
%    which skips values xmin with finite size bias > 0.1.

vec     = [];
sample  = [];
xminx   = [];
limit   = [];
finite  = false;
nosmall = false;
nowarn  = false;

% parse command-line parameters; trap for bad input
i=1; 
while i<=length(varargin) 
  argok = 1; 
  if ischar(varargin{i}) 
    switch varargin{i}
        case 'range',        vec     = varargin{i+1}; i = i + 1;
        case 'sample',       sample  = varargin{i+1}; i = i + 1;
        case 'limit',        limit   = varargin{i+1}; i = i + 1;
        case 'xmin',         xminx   = varargin{i+1}; i = i + 1;
        case 'finite',       finite  = true;
        case 'nowarn',       nowarn  = true;
        case 'nosmall',      nosmall = true;
        otherwise, argok=0; 
    end
  end
  if ~argok 
    disp(['(PLFIT) Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end
if ~isempty(vec) && (~isvector(vec) || min(vec)<=1)
	fprintf('(PLFIT) Error: ''range'' argument must contain a vector; using default.\n');
    vec = [];
end
if ~isempty(sample) && (~isscalar(sample) || sample<2)
	fprintf('(PLFIT) Error: ''sample'' argument must be a positive integer > 1; using default.\n');
    sample = [];
end
if ~isempty(limit) && (~isscalar(limit) || limit<min(x))
	fprintf('(PLFIT) Error: ''limit'' argument must be a positive value >= 1; using default.\n');
    limit = [];
end
if ~isempty(xminx) && (~isscalar(xminx) || xminx>=max(x))
	fprintf('(PLFIT) Error: ''xmin'' argument must be a positive value < max(x); using default behavior.\n');
    xminx = [];
end

% reshape input vector
x = reshape(x,numel(x),1);

% select method (discrete or continuous) for fitting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x)
    f_dattype = 'REAL';
else                 
    f_dattype = 'UNKN';
end
if strcmp(f_dattype,'INTS') && min(x) > 1000 && length(x)>100
    f_dattype = 'REAL';
end

% estimate xmin and alpha, accordingly
switch f_dattype
    
    case 'REAL'
        xmins = unique(x);
        xmins = xmins(1:end-1);
        if ~isempty(xminx)
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end
        if ~isempty(limit)
            xmins(xmins>limit) = [];
        end
        if ~isempty(sample)
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end
        dat   = zeros(size(xmins));
        z     = sort(x);
        for xm=1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin); 
            n    = length(z);
            % estimate alpha using direct MLE
            a    = n ./ sum( log(z./xmin) );
            if nosmall
                if (a-1)/sqrt(n) > 0.1
                    dat(xm:end) = [];
                    xm = length(xmins)+1; %#ok<FXSET,NASGU>
                    break;
                end
            end
            % compute KS statistic
            cx   = (0:n-1)'./n;
            cf   = 1-(xmin./z).^a;
            dat(xm) = max( abs(cf-cx) );
        end
        D     = min(dat);
        xmin  = xmins(find(dat<=D,1,'first'));
        z     = x(x>=xmin);
        n     = length(z); 
        alpha = 1 + n ./ sum( log(z./xmin) );
        if finite, alpha = alpha*(n-1)/n+1/n; end % finite-size correction
        if n < 50 && ~finite && ~nowarn
%             fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end
        L = n*log((alpha-1)/xmin) - alpha.*sum(log(z./xmin));

    case 'INTS'
        
        if isempty(vec)
            vec  = (1.50:0.01:3.50);    % covers range of most practical 
        end                            % scaling parameters
        zvec = zeta(vec);

        xmins = unique(x);
        xmins = xmins(1:end-1);
        if ~isempty(xminx)
            xmins = xmins(find(xmins>=xminx,1,'first'));
        end
        if ~isempty(limit)
            limit = round(limit);
            xmins(xmins>limit) = [];
        end
        if ~isempty(sample)
            xmins = xmins(unique(round(linspace(1,length(xmins),sample))));
        end
        if isempty(xmins)
            fprintf('(PLFIT) Error: x must contain at least two unique values.\n');
            alpha = NaN; xmin = x(1); D = NaN; %#ok<NASGU>
            return;
        end
        xmax   = max(x);
        dat    = zeros(length(xmins),2);
        z      = x;
        fcatch = 0;

        for xm=1:length(xmins)
            xmin = xmins(xm);
            z    = z(z>=xmin);
            n    = length(z);
            % estimate alpha via direct maximization of likelihood function
            if fcatch==0
                try
                    % vectorized version of numerical calculation
                    zdiff = sum( repmat((1:xmin-1)',1,length(vec)).^-repmat(vec,xmin-1,1) ,1);
                    L = -vec.*sum(log(z)) - n.*log(zvec - zdiff);
                catch
                    % catch: force loop to default to iterative version for
                    % remainder of the search
                    fcatch = 1;
                end
            end
            if fcatch==1
                % force iterative calculation (more memory efficient, but 
                % can be slower)
                L       = -Inf*ones(size(vec));
                slogz   = sum(log(z));
                xminvec = (1:xmin-1);
                for k=1:length(vec)
                    L(k) = -vec(k)*slogz - n*log(zvec(k) - sum(xminvec.^-vec(k)));
                end
            end
            [Y,I] = max(L); %#ok<ASGLU>
            % compute KS statistic
            fit = cumsum((((xmin:xmax).^-vec(I)))./ (zvec(I) - sum((1:xmin-1).^-vec(I))));
            cdi = cumsum(histcounts(discretize(z,length(xmin:xmax)))./n);
            % cdi = cumsum(hist(z,xmin:xmax)./n);
            dat(xm,:) = [max(abs( fit - cdi )) vec(I)];
        end
        % select the index for the minimum value of D
        [D,I] = min(dat(:,1)); %#ok<ASGLU>
        xmin  = xmins(I);
        z     = x(x>=xmin);
        n     = length(z);
        alpha = dat(I,2);
        if finite, alpha = alpha*(n-1)/n+1/n; end % finite-size correction
        if n < 50 && ~finite && ~nowarn
%             fprintf('(PLFIT) Warning: finite-size bias may be present.\n');
        end
        L     = -alpha*sum(log(z)) - n*log(zvec(find(vec<=alpha,1,'last')) - sum((1:xmin-1).^-alpha));

    otherwise
        fprintf('(PLFIT) Error: x must contain only reals or only integers.\n');
        alpha = [];
        xmin  = [];
        L     = [];
        return;
end