function [vslist] = i_checkexistingembed(sce, ndim)
%I_CHECKEXISTINGEMBED  Names of the embeddings SCE really carries.
%
%   Delegates to PKG.E_HASEMBEDDING so that this test has one
%   implementation. +cli and @SingleCellExperiment need the same predicate
%   and cannot reach into +gui for it; both had grown their own version
%   based on isempty(sce.s), which is never true because the constructor
%   fills S with randn.
if nargin < 2, ndim = []; end

[tf, names] = pkg.e_hasembedding(sce, ndim);
if tf
    vslist = names;
else
    vslist = '';
end
