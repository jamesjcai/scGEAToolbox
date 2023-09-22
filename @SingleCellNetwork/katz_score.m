function s = katz_score(obj, b)
%Katz score is a potential link between nodes i and j i
%ref: http://www.cs.sandia.gov/~dmdunla/publications/DuKoAc10.pdf

if nargin < 2, b = 0.8; end
A = obj.A;
n = size(A, 1);
if issparse(A)
    I = speye(n);
    s = inv(I-b*A) - I;
else
    I = eye(n);
    s = pinv(I-b*A) - I;
end
end
