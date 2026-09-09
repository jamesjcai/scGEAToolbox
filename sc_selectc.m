function [X, keptidx] = sc_selectc(X, libszcutoff, gnnumcutoff)
%SC_SELECTC Select cells by library size and number of genes
%   [Y, keptidx] = sc_selectc(X, libszcutoff, gnnumcutoff) keeps the
%   columns of the genes-by-cells matrix X whose library size and
%   detected-gene count both clear their cutoffs.
%
%   LIBSZCUTOFF has two meanings, distinguished by magnitude:
%
%     >= 1        an absolute library size, in counts (default 1000)
%     in (0, 1)   a quantile of the observed library sizes, so 0.15
%                 drops the smallest 15% of cells
%
%   GNNUMCUTOFF is always an absolute number of detected genes
%   (default 500).
%
%   The split used to be at "> 1.0", which put 1 on the quantile side.
%   quantile(libsz, 1) is the largest library in the matrix, so asking
%   for a cutoff of 1 -- the obvious way to write "cells with at least
%   one count" -- kept the single biggest cell and discarded everything
%   else, while 1.0000001 kept all of them. On a 400-cell test matrix
%   that was 1 cell against 400, and through SC_QCFILTER, whose caller
%   in GUI.CALLBACK_SELECTCELLSBYQC takes the cutoff from a text box and
%   only checks it is positive, a 200x400 matrix came back as 97x1 with
%   nothing said.
%
%   See also SC_QCFILTER, SC_FILTERC, SC_SELECTG.

arguments
    X {mustBeNumeric}
    libszcutoff (1, 1) double {mustBeNonnegative, mustBeFinite} = 1000
    gnnumcutoff (1, 1) double {mustBeNonnegative, mustBeFinite} = 500
end

libsz = full(sum(X, 1));
gnnum = full(sum(X > 0, 1));

if libszcutoff >= 1
    keptlib = libsz >= libszcutoff;
else
    keptlib = libsz >= quantile(libsz, libszcutoff);
end
keptidx = keptlib & (gnnum >= gnnumcutoff);

if ~any(keptidx) && ~isempty(libsz)
    % Returning an empty matrix is the honest answer, but silently is
    % not: SC_QCFILTER loops on the size until it stops changing, so an
    % over-strict cutoff there empties the genes as well and the caller
    % is left holding 0x0 with no clue which cutoff did it.
    warning('sc_selectc:noCellsKept', ...
        ['No cell clears both cutoffs, so no cells are kept. The ', ...
        'library-size cutoff is %g against a largest library of %g, ', ...
        'and the detected-gene cutoff is %g against a maximum of %g ', ...
        'genes detected in any one cell.'], ...
        libszcutoff, max(libsz), gnnumcutoff, max(gnnum));
end

X = X(:, keptidx);

end
