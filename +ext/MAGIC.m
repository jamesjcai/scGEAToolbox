function X=MAGIC(X,donorm)

if nargin<2, donorm=true; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','MAGIC');
addpath(pth);

% gene_names=cellstr(gl123);  MAGIC needs [cells x genes]
% data=X';

% library size normalization
% libsize = sum(data,2);
% data = bsxfun(@rdivide, data, libsize) * median(libsize);
if donorm
    X=sc_norm(X,'type','libsize');
    % log transform -- usually one would log transform the data. Here we don't do it.
    % data = log(data + 0.1);
end

% MAGIC needs [cells x genes]
X=X';

[pc_imputed, U, ~] = run_magic_original(X, 'npca', 100, 'k', 15, 'a', 15, 'make_plot_opt_t', false);

% plot_genes = {'Cdh1', 'Vim', 'Fn1', 'Zeb1'};
% [M_imputed, genes_found] = project_genes(plot_genes, gene_names, pc_imputed, U);
M = pc_imputed * U';    % project
X=M';

end