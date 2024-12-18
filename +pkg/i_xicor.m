function res = i_xicor(x, y, pvalue, ties, method, nperm)
% This function computes the xi coefficient between two vectors x and y. If only one coefficient is computed it can be used to test independence using a Monte Carlo
% permutation test or through an asymptotic approximation test.
% x Vector of numeric values in the first coordinate.
% y Vector of numeric values in the second coordinate.
% pvalue Whether or not to return the p-value of rejecting independence, if TRUE the
% function also returns the standard deviation of xi.
% ties Do we need to handle ties? If ties=TRUE the algorithm assumes that the data
% has ties and employs the more elaborated theory for calculating s.d. and P-value.
% Otherwise, it uses the simpler theory. There is no harm in putting ties = TRUE
% even if there are no ties.
% method If method = "asymptotic" the function returns P-values computed by the asymptotic theory. If method = "permutation", a permutation test with nperm permutations is employed to estimate the P-value. Usually, there is no need for the
% permutation test. The asymptotic theory is good enough.
% nperm In the case of a permutation test, nperm is the number of permutations to do.
% In the case pvalue=FALSE, function returns the value of the xi coefficient, if the input is a matrix, a
% matrix of coefficients is returned. In the case pvalue=TRUE is chosen, the function returns a list:
% xi The value of the xi coefficient.
% sd The standard deviation.
% pval The test p-value.
% References
% Chatterjee, S. (2020) <arXiv:1909.10140>.
if nargin < 3
    pvalue = false;
end
if nargin < 4
    ties = true;
end
if nargin < 5
    method = 'asymptotic';
end
if nargin < 6
    nperm = 1000;
end
x = reshape(x,[],1);
y = reshape(y,[],1);
result = calculateXI(x, y, false);
xi = result.xi;
CU = result.CU;
n = length(x);
if pvalue
    if ~ties
        res = struct('xi', xi, 'sd', sqrt(2/(5 * n)), 'pval', 1 - normcdf(sqrt(n) * xi/sqrt(2/5)));
    else
        if strcmp(method, 'asymptotic')
            fr = result.fr;
            qfr = sort(fr);
            ind = (1:n)';
            ind2 = 2 * n - 2 * ind + 1;
            ai = mean(ind2 .* qfr .* qfr)/n;
            ci = mean(ind2 .* qfr)/n;
            cq = cumsum(qfr);
            m = (cq + (n - ind) .* qfr)/n;
            b = mean(m.^2);
            v = (ai - 2 * b + ci^2)/(CU^2);
            res = struct('xi', xi, 'sd', sqrt(v/n), 'pval', 1 - normcdf(sqrt(n) * xi/sqrt(v)));
        elseif strcmp(method, 'permutation')
            rp = zeros(1, nperm);
            for i = 1:nperm
                x1 = rand(n, 1);
                rp(i) = calculateXI(x1, y, true);
            end
            res = struct('xi', xi, 'sd', sqrt(var(rp)), 'pval', mean(rp > xi));
        else
            error('method for test can only be asymptotic or permutation');
        end
    end
else
    res = xi;
end
end
function res = calculateXI(x, y, simple)
% Arguments
% x Vector of numeric values in the first coordinate.
% y Vector of numeric values in the second coordinate.
% simple Whether auxiliary information is kept to pass on.
% Value
% In the case simple = TRUE, function returns the value of the xi coefficient, If simple = FALSE is
% chosen, the function returns a list:
% xi The xi coefficient
% fr rearranged rank of yvec
% CU mean(gr*(1-gr))
% Note
% Auxiliary function with no checks for NA, etc.
% Author(s)
% Sourav Chatterjee, Susan Holmes
% References
% Chatterjee, S. (2020) A New Coefficient Of Correlation, <arXiv:1909.10140>
if nargin < 3
    simple = true;
end
n = length(x);
PI = rankR(x, 'random');
fr = rankR(y, 'max') / n;
gr = rankR(-y, 'max') / n;
[~, ord] = sort(PI);
fr = fr(ord);
A1 = sum(abs(fr(1:(n - 1)) - fr(2:n))) / (2 * n);
CU = mean(gr .* (1 - gr));
xi = 1 - A1 / CU;
if simple
    res = xi;
else
    res = struct('xi', xi, 'fr', fr, 'CU', CU);
end
end
function Rx = rankR(x, ties_method)
n = numel(x);
Rx = zeros(size(x));
switch ties_method
    case 'random'
        xx = x + rand(size(x))*1e-10;
        [~ ,r] = sort(xx);
        Rx(r) = 1:n;
    case 'max'
        [~ ,r] = sort(x);
        Rx(r) = 1:n;
        for i = 1:n
            position = ( x == x(i) );
            Rx(position) = max(Rx(position));
        end
end
end
