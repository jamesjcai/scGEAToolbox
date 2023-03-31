function [eigvector, eigvalue] = KPCA(data, options)
%KPCA	Kernel Principal Component Analysis
%
%	Usage:
%       [eigvector, eigvalue] = KPCA(data, options)
% 
%             Input:
%               data    - 
%                      if options.Kernel = 0
%                           Data matrix. Each row vector of fea is a data
%                           point. 
%                      if options.Kernel = 1
%                           Kernel matrix. 
%
%             options   - Struct value in Matlab. The fields in options
%                         that can be set:
%                      Kernel  -  1: data is actually the kernel matrix. 
%                                 0: ordinary data matrix. 
%                                   Default: 0 
%                         
%                       Please see constructKernel.m for other Kernel options. 
%
%                     ReducedDim   - The dimensionality of the reduced subspace. If 0,
%                         all the dimensions will be kept. Default is 30. 
%
%             Output:
%               eigvector - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = K(x,:)*eigvector
%                           will be the embedding result of x.
%                           K(x,:) = [K(x1,x),K(x2,x),...K(xm,x)]
%               eigvalue  - The sorted eigvalue of PCA eigen-problem. 
%
%	Examples:
%           options.KernelType = 'Gaussian';
%           options.t = 1;
%           options.ReducedDim = 4;
% 			fea = rand(7,10);
% 			[eigvector,eigvalue] = KPCA(fea,options);
%           feaTest = rand(3,10);
%           Ktest = constructKernel(feaTest,fea,options)
%           Y = Ktest*eigvector;
% 
%Reference:
%
%   Bernhard Sch�lkopf, Alexander Smola, Klaus-Robert M�ller, �Nonlinear
%   Component Analysis as a Kernel Eigenvalue Problem", Neural Computation,
%   10:1299-1319, 1998.
%
% 
%   version 1.1 --Dec./2011 
%   version 1.0 --April/2005 
%
%   Written by Deng Cai (dengcai AT gmail.com)
%                                                   
MAX_MATRIX_SIZE = 1600; % You can change this number according your machine computational power
EIGVECTOR_RATIO = 0.1; % You can change this number according your machine computational power


ReducedDim = 30;
if isfield(options,'ReducedDim')
    ReducedDim = options.ReducedDim;
end


if isfield(options,'Kernel') && options.Kernel
    K = data;
else
    K = constructKernel(data,[],options);
end
clear data;
    
nSmp = size(K,1);
if (ReducedDim > nSmp) || (ReducedDim <=0)
    ReducedDim = nSmp;
end

sumK = sum(K,2);
H = repmat(sumK./nSmp,1,nSmp);
K = K - H - H' + sum(sumK)/(nSmp^2);
K = max(K,K');
clear H;

if nSmp > MAX_MATRIX_SIZE && ReducedDim < nSmp*EIGVECTOR_RATIO
    % using eigs to speed up!
    option = struct('disp',0);
    [eigvector, eigvalue] = eigs(K,ReducedDim,'la',option);
    eigvalue = diag(eigvalue);
else
    [eigvector, eigvalue] = eig(K);
    eigvalue = diag(eigvalue);
    
    [~, index] = sort(-eigvalue);
    eigvalue = eigvalue(index);
    eigvector = eigvector(:,index);
end


if ReducedDim < length(eigvalue)
    eigvalue = eigvalue(1:ReducedDim);
    eigvector = eigvector(:, 1:ReducedDim);
end

maxEigValue = max(abs(eigvalue));
eigIdx = find(abs(eigvalue)/maxEigValue < 1e-6);
eigvalue (eigIdx) = [];
eigvector (:,eigIdx) = [];

for i=1:length(eigvalue) % normalizing eigenvector
    eigvector(:,i)=eigvector(:,i)/sqrt(eigvalue(i));
end


end

function K = constructKernel(fea_a,fea_b,options)
% function K = constructKernel(fea_a,fea_b,options)
%	Usage:
%	K = constructKernel(fea_a,[],options)
%
%   K = constructKernel(fea_a,fea_b,options)
%
%	fea_a, fea_b  : Rows of vectors of data points. 
%
%   options       : Struct value in Matlab. The fields in options that can
%                   be set: 
%           KernelType  -  Choices are:
%               'Gaussian'      - e^{-(|x-y|^2)/2t^2}
%               'Polynomial'    - (x'*y)^d
%               'PolyPlus'      - (x'*y+1)^d
%               'Linear'        -  x'*y
%
%               t       -  parameter for Gaussian
%               d       -  parameter for Poly
%
%   version 1.0 --Sep/2006 
%
%   Written by Deng Cai (dengcai2 AT cs.uiuc.edu)
%

if (~exist('options','var'))
   options = [];
else
   if ~isstruct(options) 
       error('parameter error!');
   end
end



%=================================================
if ~isfield(options,'KernelType')
    options.KernelType = 'Gaussian';
end

switch lower(options.KernelType)
    case {lower('Gaussian')}        %  e^{-(|x-y|^2)/2t^2}
        if ~isfield(options,'t')
            options.t = 1;
        end
    case {lower('Polynomial')}      % (x'*y)^d
        if ~isfield(options,'d')
            options.d = 2;
        end
    case {lower('PolyPlus')}      % (x'*y+1)^d
        if ~isfield(options,'d')
            options.d = 2;
        end
    case {lower('Linear')}      % x'*y
    otherwise
        error('KernelType does not exist!');
end


%=================================================

switch lower(options.KernelType)
    case {lower('Gaussian')}       
        if isempty(fea_b)
            D = EuDist2(fea_a,[],0);
        else
            D = EuDist2(fea_a,fea_b,0);
        end
        K = exp(-D/(2*options.t^2));
    case {lower('Polynomial')}     
        if isempty(fea_b)
            D = full(fea_a * fea_a');
        else
            D = full(fea_a * fea_b');
        end
        K = D.^options.d;
    case {lower('PolyPlus')}     
        if isempty(fea_b)
            D = full(fea_a * fea_a');
        else
            D = full(fea_a * fea_b');
        end
        K = (D+1).^options.d;
    case {lower('Linear')}     
        if isempty(fea_b)
            K = full(fea_a * fea_a');
        else
            K = full(fea_a * fea_b');
        end
    otherwise
        error('KernelType does not exist!');
end

if isempty(fea_b)
    K = max(K,K');
end    
end

function D = EuDist2(fea_a,fea_b,bSqrt)
%EUDIST2 Efficiently Compute the Euclidean Distance Matrix by Exploring the
%Matlab matrix operations.
%
%   D = EuDist(fea_a,fea_b)
%   fea_a:    nSample_a * nFeature
%   fea_b:    nSample_b * nFeature
%   D:      nSample_a * nSample_a
%       or  nSample_a * nSample_b
%
%    Examples:
%
%       a = rand(500,10);
%       b = rand(1000,10);
%
%       A = EuDist2(a); % A: 500*500
%       D = EuDist2(a,b); % D: 500*1000
%
%   version 2.1 --November/2011
%   version 2.0 --May/2009
%   version 1.0 --November/2005
%
%   Written by Deng Cai (dengcai AT gmail.com)


if ~exist('bSqrt','var')
    bSqrt = 1;
end

if (~exist('fea_b','var')) || isempty(fea_b)
    aa = sum(fea_a.*fea_a,2);
    ab = fea_a*fea_a';
    
    if issparse(aa)
        aa = full(aa);
    end
    
    D = bsxfun(@plus,aa,aa') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
    D = max(D,D');
else
    aa = sum(fea_a.*fea_a,2);
    bb = sum(fea_b.*fea_b,2);
    ab = fea_a*fea_b';

    if issparse(aa)
        aa = full(aa);
        bb = full(bb);
    end

    D = bsxfun(@plus,aa,bb') - 2*ab;
    D(D<0) = 0;
    if bSqrt
        D = sqrt(D);
    end
end

end