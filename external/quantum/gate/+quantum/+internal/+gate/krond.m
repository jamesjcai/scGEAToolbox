function K = krond(varargin)
% Same as KRON, but takes any number of inputs

%   Copyright 2021-2022 The MathWorks, Inc.

arguments(Repeating)
    varargin {mustBeNonsparse} 
end

d = length(varargin);

for ii=1:d
    if ~ismatrix(varargin{ii})
        error(message('MATLAB:kron:TwoDInput'));
    end
end

if nargin <= 5
    if nargin == 0
        K = [];
    elseif nargin == 1
        K = varargin{1};
    elseif nargin == 2
        K = locKron(varargin{1}, varargin{2});
    elseif nargin == 3
        K = locKron(varargin{1}, varargin{2});
        K = locKron(K, varargin{3});
    elseif nargin == 4
        K = locKron(varargin{1}, varargin{2});
        K3 = locKron(varargin{3}, varargin{4});
        K = locKron(K, K3);
    else % nargin == 5
        K = locKron(varargin{1}, varargin{2});
        K3 = locKron(varargin{3}, varargin{4});
        K = locKron(K, K3);
        K = locKron(K, varargin{5});
    end
    return
end

% Binary tree approach: First multiply together every pair of matrices,
% then multiply those results (pairwise again) until only the first matrix
% remains, a combination of all operations.
p = 1;
while p < d
    for ii=1:2*p:d-p
        varargin{ii} = locKron(varargin{ii}, varargin{ii+p});
        varargin{ii+p} = [];
    end
    p = 2*p;
end

K = varargin{1};

% Direct approach (similar performance to just calling KRON in a loop):
% K = 1;
% for ii=1:d
%     K = locKron(K, varargin{ii});
% end

end

function K = locKron(K, A)
    mK = size(K, 1);
    nK = size(K, 2);
    mA = size(A, 1);
    nA = size(A, 2);
    
    K = reshape(K, 1, mK, 1, nK);
    A = reshape(A, mA, 1, nA, 1);
    K = K.*A;
    
    K = reshape(K, mK*mA, nK*nA);
end
