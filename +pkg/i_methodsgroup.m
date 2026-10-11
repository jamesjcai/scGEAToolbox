function group = i_methodsgroup(names, step)
%I_METHODSGROUP Stored Methods-text keys that belong to one step.
%
%   group = pkg.i_methodsgroup(names, step) returns the entries of NAMES
%   equal to STEP or starting with STEP + "_", e.g. "embedding_features"
%   and "embedding_umap3d" for "embedding". A "_features" sub-key comes
%   first, because it describes what the other sub-keys were computed on;
%   the rest keep their stored order.
%
%   See also SC_METHODSTEXT, PKG.I_METHODSSHOW.

names = string(names(:))';
group = names(names == step | startsWith(names, step + "_"));
isfeatures = endsWith(group, "_features");
group = [group(isfeatures), group(~isfeatures)];
end
