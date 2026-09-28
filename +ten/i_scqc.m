function [X, genelist] = i_scqc(X, genelist, opts)
%I_SCQC  Quality control and CPM normalization, as the R scTenifold packages.
%   [X, GENELIST] = TEN.I_SCQC(X, GENELIST) applies scQC() and then
%   cpmNormalization() of scTenifoldNet/scTenifoldKnk to the genes-by-cells
%   count matrix X, and returns the normalized matrix with the genes that
%   survive, in their input order. The steps, all with R's strict
%   inequalities:
%     1. negative values set to 0
%     2. cells with library size > MinLibSize
%     3. cells whose library size is not an outlier of boxplot.stats (1.5
%        interquartile ranges beyond Tukey's hinges)
%     4. cells with mitochondrial (^MT-) ratio < MaxMTRatio, if there are
%        mitochondrial genes
%     5. genes expressed in more than MinPCT of the remaining cells
%     6. CPM: each cell divided by its library size, times 1e6 (no log)
%
%   Used by the reference mode of TEN.SCTENIFOLDNET and TEN.SCTENIFOLDKNK.
%
% see also: TEN.SCTENIFOLDNET, TEN.SCTENIFOLDKNK

arguments
    X {mustBeNumeric}
    genelist string
    opts.MinLibSize (1, 1) double = 1000
    opts.RemoveOutlierCells (1, 1) logical = true
    opts.MinPCT (1, 1) double = 0.05
    opts.MaxMTRatio (1, 1) double = 0.1
end

X = full(double(X));
genelist = genelist(:);
X(X < 0) = 0;

X = X(:, sum(X, 1) > opts.MinLibSize);

if opts.RemoveOutlierCells && ~isempty(X)
    libSize = sum(X, 1);
    hinges = i_fivenum(libSize);
    iqrWidth = hinges(4) - hinges(2);
    isOutlier = libSize < hinges(2) - 1.5*iqrWidth | libSize > hinges(4) + 1.5*iqrWidth;
    % R drops every cell whose library size equals an outlier value
    X = X(:, ~ismember(libSize, libSize(isOutlier)));
end

isMito = startsWith(upper(genelist), "MT-");
if any(isMito)
    mtRatio = sum(X(isMito, :), 1)./sum(X, 1);
    X = X(:, mtRatio < opts.MaxMTRatio);
end

keepGene = mean(X ~= 0, 2) > opts.MinPCT;
X = X(keepGene, :);
genelist = genelist(keepGene);

% Same order of operations as cpmNormalization
X = X./sum(X, 1)*1e6;
end

function h = i_fivenum(x)
% Tukey's five-number summary, as fivenum() in R
x = sort(x(:));
n = numel(x);
n4 = floor((n + 3)/2)/2;
d = [1, n4, (n + 1)/2, n + 1 - n4, n];
h = 0.5*(x(floor(d)) + x(ceil(d)));
end
