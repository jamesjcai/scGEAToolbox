function g = i_defaultgenenames(n, prefix)
%I_DEFAULTGENENAMES  Generate a default gene-name string column vector.
%
%   g = pkg.i_defaultgenenames(n)         returns ["Gene1"; "Gene2"; ...]
%   g = pkg.i_defaultgenenames(n, prefix) uses the given prefix string
%
%   Example:
%       g = pkg.i_defaultgenenames(5);          % "Gene1".."Gene5"
%       g = pkg.i_defaultgenenames(5, "G");     % "G1".."G5"

if nargin < 2, prefix = "Gene"; end
g = prefix + string((1:n).');
end
