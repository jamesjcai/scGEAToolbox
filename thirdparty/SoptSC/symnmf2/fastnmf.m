function [A,B] = fastnmf(X,N,varargin)
% FASTNMF Fast least-squares non-negative matrix factorization 
%  Approximates X by A*B s.t. A,B>=0
%
% Usage
%   [A,B] = fastnmf(X,N,[options])
%
% Input
%   X           Data matrix (I x J)
%   N           Number of components 
%   options
%     maxiter   Number of iterations (default 100)
%     iiter     Number of inner iterations (default 10)
%     A         Initial value for A (I x N)
%     B         Initial value for B (N x J)
%
% Output
%   A           Samples of A (I x N x M)
%   B           Samples of B (N x J x M)
%
% Author
%   Mikkel N. Schmidt, 
%   DTU Informatics, Technical University of Denmark. 

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

[I,J] = size(X);

opts = mgetopt(varargin);
maxiter = mgetopt(opts, 'maxiter', 100);
iiter = mgetopt(opts, 'iiter', 10);
A = mgetopt(opts, 'A', rand(I,N));
B = mgetopt(opts, 'B', rand(N,J));

tol = 1e-8;

for m = 1:maxiter
    BBt = B*B';
    XBt = X*B';
    cA = true;
    r = 0;
	while cA & r<iiter
        r = r+1;
        gA = A*BBt-XBt;    
        aA = A~=0 | gA<0;
        cA = find(~all(abs(min(gA,A))<tol,2));
        [uc, jc, ic] = unique(aA(cA,:), 'rows');
        for j = 1:length(jc)
            c = uc(j,:);
            i = cA(ic==j);
            A(i,c) = max(XBt(i,c)/BBt(c,c),0);
        end
    end
     
    AtA = A'*A;
    AtX = A'*X;
    cB = true;
    r = 0;
	while cB & r<iiter
        r = r+1;
        gB = AtA*B-AtX;
        aB = B~=0 | gB<0;
        cB = find(~all(abs(min(gB,B))<tol,1));
        [uc, ic, jc] = unique(aB(:,cB)','rows');
        for i = 1:length(ic)
            c = uc(i,:)';
            j = cB(jc==i);
            B(c,j) = max(AtA(c,c)\AtX(c,j),0);
        end  
    end

%     KKT = max(max([abs(min((A*B-X)*B',A)); abs(min(A'*(A*B-X),B))']))
%     norm(X-A*B,'fro')
end

%--------------------------------------------------------------------------
function out = mgetopt(varargin)
% MGETOPT Parser for optional arguments
% 
% Usage
%   Get alpha parameter structure from 'varargin'
%     opts = mgetopt(varargin);
%
%   Get and parse alpha parameter:
%     var = mgetopt(opts, varname, default);
%        opts:    parameter structure
%        varname: name of variable
%        default: default value if variable is not set
%
%     var = mgetopt(opts, varname, default, command, argument);
%        command, argument:
%          String in set:
%          'instrset', {'str1', 'str2', ... }
%
% Example
%    function y = myfun(x, varargin)
%    ...
%    opts = mgetopt(varargin);
%    parm1 = mgetopt(opts, 'parm1', 0)
%    ...

% Copyright 2007 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

if nargin==1
    if isempty(varargin{1})
        out = struct;
    elseif isstruct(varargin{1})
        out = varargin{1}{:};
    elseif isstruct(varargin{1}{1})
        out = varargin{1}{1};
    else
        out = cell2struct(varargin{1}(2:2:end),varargin{1}(1:2:end),2);
    end
elseif nargin>=3
    opts = varargin{1};
    varname = varargin{2};
    default = varargin{3};
    validation = varargin(4:end);
    if isfield(opts, varname)
        out = opts.(varname);
    else
        out = default;
    end
    
    for narg = 1:2:length(validation)
        cmd = validation{narg};
        arg = validation{narg+1};
        switch cmd
            case 'instrset'
                if ~any(strcmp(arg, out))
                    fprintf(['Wrong argument %sigma = ''%sigma'' - ', ...
                        'Using default : %sigma = ''%sigma''\n'], ...
                        varname, out, varname, default);
                    out = default;
                end
            case 'dim'
                if ~all(size(out)==arg)
                    fprintf(['Wrong argument dimension: %sigma - ', ...
                        'Using default.\n'], ...
                        varname);
                    out = default;
                end
            otherwise
                error('Wrong option: %sigma.', cmd);
        end
    end
end