function [value, functionHandle] = splineMI(x,y,M,k)
%SPLINEMI - a function that calculates approximate mutual
%information between an Nx1 vector, x, and an NxP matrix, y. The algorithm
%calculates mutual information using binning and splines as described in
%[1]. 
%
% The estimator in [1] is subject to bias. Choice of M and k strongly
% affect this bias in ways discussed in [1], but also, depending on the
% underlying distribution of x and y, the number of samples from x and y
% affect this bias. For the estimate of mutual information to converge in
% a predictable way, the underlying distribution of x and y should be bounded
% above and below. 
%
% It's also important to realize that the estimator converges to the mutual
% information of two discrete, proxy random variables for x and y that are
% defined by the binning process. If x and y are sampled from a continuous
% random variable, the estimated mutual information is not, in general, the
% mutual information between the two continuous random variables.
%
% Usage:
% bins = 6;
% spline_order = 3;
% x = linspace(-1,1);
% y = x.^2 + .1*randn(size(x))
% mi = splinemi(x,y,bins,spline_order);
%
% If y has many columns, the most expensive component of this calculation
% is calculating the entropy and binning of the variable y. For many
% applications, it is necessary to calculate the mutual information between
% several different x values and the same y. To accelerate this use case we
% provide a function handle that takes new values of x, and returns the
% mutual information between the new x and the original y. Example:
%
% bins = 6;
% spline_order = 3;
% x = linspace(-1,1,100)';
% noisyMatrix = randn(100,100);
% noisyMatrix = bsxfun(@times,linspace(0,1),noisyMatrix);
% y = bsxfun(@plus,x.^2,noisyMatrix);
% [mi, func] = splinemi(x,y,bins,splin_order);
%
% for iter = 1:10
%    miDegree(iter) = func(x.^iter); % gives the mutual information between y and x.^iter
% end
%
% note: mi is equal to miDegree(2) in the example above.
%
% References: 
% [1] Daub, Carsten O., et al. "Estimating mutual information using
%     B-spline functions--an improved similarity measure for analysing gene
%     expression data." BMC bioinformatics 5.1 (2004): 118.

% Copyright MathWorks 2014

%% Validate inputs

validateY(y);
[numRows, numCols] = size(y);
if ~isempty(x)
    validateX(x,numRows);
end

if ~(isscalar(M) && floor(M) == M)
    error('bioinfo:mutualinformation:Minteger','M, the number of bins, must be an integer');
end
if ~(isscalar(k) &&floor(k) == k)
    error('bioinfo:mutualinformation:Kinteger','k, the spline order, must be an integer');
end
if M<k
    error('bioinfo:mutualinformation:BinsSpline','The number of bins, M, must be greater than or equal to k, the spline order');
end

%% get knot as function of k/M

t = [zeros(1,k-1), 0:(M-k+1) , (M-k+1)*ones(1,k-1)];

%% Calculate entropy

entropy_y = zeros(1,numCols);
bval_y    = zeros(numRows,M,numCols);
for ind = 1:numCols
    [entropy_y(ind), bval_y(:,:,ind)] = bioinfo.internal.entropy(y(:,ind),k,t);
end

bval_y_stacked = reshape(bval_y,[numRows, numCols*M]);

%% calculate joint prob/entropy

small = eps(1e-7);
[~, numCols] = size(y);

functionHandle = @calculateMutualInfoFromEntropy;

if isempty(x)
    value = [];
else
    value = functionHandle(x);
end
    
    function mutual_info = calculateMutualInfoFromEntropy(x)
        vecLength = length(x);
        
        [entropy_x, bval_x] = bioinfo.internal.entropy(x,k,t);

        joint_prob = bval_x'*bval_y_stacked/vecLength; % perform one matrix multiple for all indices, then index into the results O(bkn)
        entropy_components = sum(joint_prob.*log(joint_prob+small)); % add small to avoid having to nansum, does not introduce error when prob <1e-7, small otherwise
        
        entropy_joint = zeros(1,numCols);
        
        for start = 1:M
            % sum indices 1:M, (M+1):2M, etc in an efficient way
            entropy_joint = entropy_joint - entropy_components(start:M:end);
        end
        
        mutual_info = entropy_x + entropy_y - entropy_joint;
        
    end

end

function bool = validateX(x,p)
bool = true;
if ~(isnumeric(x) && isreal(x) && iscolumn(x))
    error('bioinfo:mutualInformation:InvalidX','x must be a real column vector');
end

if ~(size(x,1) == p)
    error('bioinfo:mutualInformation:InvalidXSize','x must be a column vector of size with the same number of rows as y.')
end

end

function bool = validateY(y)
    bool = true;
    if ~(isnumeric(y) && isreal(y) && ismatrix(y))
        error('bioinfo:mutualInformation:InvalidY','Y must be a real matrix.')
    end
end

