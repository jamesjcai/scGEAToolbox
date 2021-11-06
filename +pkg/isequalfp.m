function yn = isequalfp(a,b)
% ISEQUALFP  Check two values for equality within floating point precision
%
% It is widely known that floating point computation has a fundamental
%   limitation: not every value can be represented exactly.  This can
%   lead to surprising results for those unfamiliar with this
%   limitation, especially since 'double' is MATLAB's default numerical
%   data type.
%
% This function accepts two float values (single or double) or arrays
%   of floats, and returns a logical value indicating whether they
%   are equal within floating point precision.  Mixed single and double
%   inputs will be evaluated based on single floating point precision.
%
% Floating point accuracy reference:
%   http://blogs.mathworks.com/loren/2006/08/23/a-glimpse-into-floating-point-accuracy/
%
% Usage:
%   yn = isequalfp(a,b)
%
%     a,b: floats or arrays of floats to compare
%
%      yn: logical scalar result indicating equality
%
% Example:
%   a = 0.3;
%   b = 0.1*3;
%   isequal(a,b)     % ans = 0
%   isequalfp(a,b)   % ans = 1
%   c = a+2*eps(a)   % c = 0.3000...
%   isequalfp(a,c)   % ans = 0
%
% See also: EPS, ISEQUAL

% v0.2 (May 2012) by Andrew Davis (addavis@gmail.com)
% todo: isequalfp(a,b,c,...)

%% Check arguments
narginchk(2,2);
assert(isfloat(a) && isfloat(b), 'inputs a and b must be floats');
assert(all(size(a) == size(b)), 'inputs a and b must be the same size');

a = a(:);   % column vectors
b = b(:);

%% Check for equivalence of each element, within the tolerance
yn = abs(a - b) <= eps(max(abs(a), abs(b)));

yn = all(yn);  % scalar logical output
