function [C] = ml_SinNLRR(X, k)
% SinNLRR -
%
% USAGE:
% >> % [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% load('example_data/example10xdata.mat');
% [C,s]=run_simlr(X,[],true);
% figure;
% scatter(s(:,1),s(:,2),20,C,'filled')

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'external', 'ml_SinNLRR');
if ~(ismcc || isdeployed)
    addpath(pth);
end
if nargin < 2 || isempty(k)
    k = fun_num_cluster(X);
    fprintf('k=%d\n', k);
end

% if nargin<3
%     donorm=false;
% end

% if donorm
%     % [X]=sc_norm(X);
%     % X=log10(X+1);
%     X=X./vecnorm(X);
% end
% The bundled external/ml_SinNLRR/SpectralClustering.m line 16 is a bare
% "warning off" with nothing to undo it, so a single call left warnings
% disabled for the rest of the MATLAB session. That is not a cosmetic
% problem: it silences every later warning the user relies on, and it made
% eleven unrelated tests in this repository's own suite fail, because every
% test that asserts a warning stopped seeing one. Restoring here rather
% than editing the third-party file keeps the fix in code we own.
warnState = warning();
restoreWarn = onCleanup(@() warning(warnState));

[C] = SinNLRRori(X, k);


end
