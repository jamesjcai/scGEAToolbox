function [xi, p] = e_xicor(x, y, varargin)
% XICOR Computes Chaterjee's xi correlation between x and y variables
%
%   [xi, p] = xicor(x, y)
%   Returns the xi-correlation with the corresponding p-value for the pair
%   of variables x and y.
%
%   Input arguments:
%
%   'x'              Independent variable. Numeric 1D array.
%
%   'y'              Dependent variable. Numeric 1D array.
%
%
%   Optional third argument:
%
%   symmetric        If true xi is computed as (r(x,y)+r(y,x))/2. Passed
%                    positionally, e.g. pkg.e_xicor(x, y, true).
%                    Default: false.
%
%   Name-value arguments:
%
%   'p_val_method'   Method to be used to compute the p-value.
%                    Options: 'theoretical' or 'permutation'.
%                    Default: 'theoretical'.
%
%   'n_perm'         Number of permutations when p_val_method is
%                    'permutation'.
%                    Default: 1000.
%
%
%   Output arguments:
%
%   'xi'             Computed xi-correlation.
%
%   'p'              Estimated p-value.
%
%
%   Notes
%   -----
%   This is an independent implementation of the method largely based on
%   the R-package developed by the original authors [3].
%   The xi-correlation is not symmetric by default.
%   Check [2] for a potential improvement over the current implementation.
%
%
%   References
%   ----------
%   [1]  Sourav Chatterjee, A New Coefficient of Correlation, Journal of
%   the American Statistical Association, 116:536, 2009-2022, 2021.
%   DOI: 10.1080/01621459.2020.1758115
%
%   [2] Zhexiao Lin* and Fang Han†, On boosting the power of Chatterjee’s
%   rank correlation, arXiv, 2021. https://arxiv.org/abs/2108.06828
%
%   [3] XICOR R package.
%   https://cran.r-project.org/web/packages/XICOR/index.html
%
%
%   Example
%   ---------
%   % Compute the xi-correlation between two variables
%
%     x = linspace(-10,10,50);
%     y = x.^2 + randn(1,50);
%     [xi, p] = xicor(x,y);
%
%
%   David Romero-Bascones, dromero@mondragon.edu
%   Biomedical Engineering Department, Mondragon Unibertsitatea, 2022
if nargin == 1
    error('err1:MoreInputsRequired', 'xicor requires at least 2 inputs.');
end
parser = inputParser;
addRequired(parser, 'x');
addRequired(parser, 'y');
% SYMMETRIC stays positional because +net/xicornet.m passes it that way,
% and the help below now says so. P_VAL_METHOD and N_PERM become real
% name-value arguments: the help has always described all three as
% name-value, but they were addOptional, so e_xicor(x, y, 'p_val_method',
% 'permutation') -- the documented call -- set SYMMETRIC to the char vector
% 'p_val_method' and died on "symmetric must be true or false".
addOptional(parser, 'symmetric', false);
addParameter(parser, 'p_val_method', 'theoretical');
addParameter(parser, 'n_perm', 1000);
parse(parser,x,y,varargin{:});
x = parser.Results.x;
y = parser.Results.y;
symmetric = parser.Results.symmetric;
p_val_method = parser.Results.p_val_method;
n_perm = parser.Results.n_perm;
if ~isnumeric(x) || ~isnumeric(y)
    error('err2:TypeError', 'x and y are must be numeric.');
end
n = length(x);
if n ~= length(y)
    error('err3:IncorrectLength', 'x and y must have the same length.');
end
if ~islogical(symmetric)
    error('err2:TypeError', 'symmetric must be true or false.');
end
% Check for NaN values
is_nan = isnan(x) | isnan(y);
if sum(is_nan) == n
    warning('No points remaining after excluding NaN.');
    xi = nan;
    return
elseif sum(is_nan) > 0
    warning('NaN values encountered.');
    x = x(~is_nan);
    y = y(~is_nan);
    n = length(x);
end
if n < 10
    warning(['Running xicor with only ', num2str(n),...
             ' points. This might produce unstable results.']);
end
[xi, r, l] = compute_xi(x, y);
if symmetric
    xi = (xi + compute_xi(y, x))/2;
end
% If only one output return xi
if nargout <= 1
    return
end
if ~strcmp(p_val_method, 'permutation') && symmetric==true
    error('err2:TypeError', ...
          'p_val_method when symmetric==true must be permutation.');
end
% Compute p-values (only valid for large n)
switch p_val_method
    case 'theoretical'
        if length(unique(y)) == n
            p = 1 - normcdf(sqrt(n)*xi, 0, sqrt(2/5));
        else
            u = sort(r);
            v = cumsum(u);
            i = 1:n;

            a = 1/n^4 * sum((2*n -2*i +1) .* u.^2);
            b = 1/n^5 * sum((v + (n - i) .* u).^2);
            c = 1/n^3 * sum((2*n -2*i +1) .* u);
            d = 1/n^3 * sum(l .* (n - l));

            tau = sqrt((a - 2*b + c^2)/d^2);

            p = 1 - normcdf(sqrt(n)*xi, 0, tau);
        end
    case 'permutation'
        xi_perm = nan(1, n_perm);

        if symmetric
            for i_perm=1:n_perm
                x_perm = x(randperm(n));
                xi_perm(i_perm) = (compute_xi(x_perm, y) + ...
                                   compute_xi(y, x_perm))/2;
            end
        else
            for i_perm=1:n_perm
                xi_perm(i_perm) = compute_xi(x(randperm(n)), y);
            end
        end

        p = sum(xi_perm > xi)/n_perm;
    otherwise
        error("Wrong p_value_method. Use 'theoretical' or 'permutation'");
end
function [xi, r, l] = compute_xi(x,y)
n = length(x);
% Reorder based on x
[~, si] = sort(x, 'ascend');
y = y(si);
% Compute y ranks
[~, si] = sort(y, 'ascend');
r = 1:n;
r(si) = r;
% If no Y ties compute it directly
if length(unique(y)) == n
    xi = 1 - 3*sum(abs(diff(r)))/(n^2 - 1);
    r = nan;
    l = nan;
else
    % Tie-corrected ranks, taken straight from the definition
    % (Chatterjee 2021): r_i = #{j : y_j <= y_i} and
    % l_i = #{j : y_j >= y_i}.
    %
    % This replaces a loop that found its tied groups with
    %
    %     idx_tie = find(groupcounts(y) > 1);
    %
    % GROUPCOUNTS on a ROW vector treats the whole row as ONE group and
    % returns a scalar 1, so find(... > 1) was empty and the loop never
    % executed -- leaving r as the plain 1:n rank and l as n - r + 1,
    % neither of which is the tie-corrected count. Every in-toolbox
    % caller passes rows: +net/xicornet.m slices X(k, :) out of a
    % genes-by-cells matrix.
    %
    % It matters most on exactly the data this toolbox handles, where
    % ties are the rule. Measured on Poisson counts with 5 distinct
    % values in 200 cells, against the definition computed directly:
    %
    %     definition   row input   column input
    %      0.005183     0.121053     0.005183
    %     -0.005571     0.113161    -0.005571
    %     -0.027100     0.162763    -0.027100
    %      0.092203     0.174230     0.092203
    %
    % The column path was already right, so the two orientations
    % disagreed by more than twentyfold. Worse, the bias has a sign:
    % for independent variables the definition scatters around zero
    % while the row path returned 0.11 to 0.17 every time, so
    % sc_grn(X, 'xicor') read a dependence into every gene pair.
    % Untied data was unaffected, which is the other branch above.
    [~, ~, grp] = unique(y(:));
    counts = accumarray(grp, 1);
    leq = cumsum(counts);              % #{y <= each distinct value}
    geq = n - [0; leq(1:end-1)];       % #{y >= each distinct value}
    r = reshape(leq(grp), 1, []);
    l = reshape(geq(grp), 1, []);

    % Compute correlation
    xi = 1 - n*sum(abs(diff(r)))/(2*sum((n - l) .* l));
end
