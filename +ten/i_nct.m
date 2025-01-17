function [XM] = i_nct(X, ptime, nsubsmpl, ncom, csubsmpl, savegrn, mksparse)
% NCT - network construction with cells subsampled using pseudotime
%
% input: X -  n (genes/features) x m (cells/samples) matrix
% input: ptime - m x 1 vector with pseudotime of cells
% output XM - k multi-layer network array (n x n x k)
import ten.*

if nargin < 7, mksparse = true; end
if nargin < 6, savegrn = true; end
if nargin < 5, csubsmpl = 500; end % number of cells in subsamples
if nargin < 4, ncom = 3; end % number of components for PC regression
if nargin < 3, nsubsmpl = 10; end % number of subsamples


[~, t] = sort(ptime);
[n, m] = size(X);
assert(length(ptime) == m)
assert(max(t) == m)
X = X(:, t);
if m < csubsmpl * 1.5, error('Too few cells.'); end

r = round(m/nsubsmpl); % r = step length
winsize = max([r, csubsmpl]);
startptx = 1:r:m;
c = 0;
while startptx(end) + winsize > m && r > 1
    c = c + 1;
    r = r - 1;
    winsize = max([r, csubsmpl]);
    startptx = 1:r:m;
    startptx = startptx(1:nsubsmpl);
end

XM = zeros(n, n, nsubsmpl);
for k = 1:nsubsmpl
    fprintf('network...%d of %d\n', k, nsubsmpl);
    Xrep = X(:, startptx(k):startptx(k)+winsize-1);
    Xrep = sc_norm(Xrep, 'type', 'libsize');
    Xrep = log1p(Xrep);
    A = sc_pcnetpar(Xrep, ncom, true);
    if savegrn
        [~, b] = fileparts(tempname);
        if ~exist(sprintf('A%d_%s.mat', k, b(1:10)), 'file')
            b = b(1:10);
        end
        save(sprintf('A%d_%s', k, b), 'A');
    end
    if mksparse
        XM(:, :, k) = ten.e_filtadjc(A, 0.95, false);
    else
        XM(:, :, k) = A;
    end
end
end
