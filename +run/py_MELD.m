function [sample_likelihoods, T] = py_MELD(X, batchid)
% MELD - a graph signal processing tool used to smooth a binary variable on
% the cell-cell graph to determine which regions of its underlying
% data manifold are enriched or depleted in cells with a specific
% feature.
arguments
    X(:, :) {mustBeNumeric}
    batchid(1, :) {mustBePositive, mustBeInteger}
end
sample_likelihoods = [];
T = [];

isdebug = false;

oldpth = pwd();
prgfoldername = 'py_MELD';

[pyok, wrkpth, x] = run.pycommon(prgfoldername);
if ~pyok, return; end
tmpfilelist = {'batchid.txt', 'input.txt', 'output.txt', 'input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
X = sc_norm(X);
X = sqrt(X)';
save('input.mat', 'X', 'batchid', '-v7.3');
disp('Input file written.');

[status] = run.pycommon2(x, wrkpth, prgfoldername);

if status == 0 && exist('output.txt', 'file')
    T = readtable('output.txt', "ReadVariableNames", true, ...
        'VariableNamingRule', 'modify');
    sample_likelihoods = table2array(T);
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
