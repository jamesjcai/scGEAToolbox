function s = ml_PHATE(X, ndim, plotit, bygene, genelist)
% RUN_PHATE
%
% PHATE is a data reduction method specifically designed for visualizing
% **high** dimensional data in **low** dimensional spaces. Details on
% this package can be found here https://github.com/KrishnaswamyLab/PHATE.
% For help, visit https://krishnaswamylab.org/get-help. For a more in
% depth discussion of the mathematics underlying PHATE, see the bioRxiv
% paper here: https://www.biorxiv.org/content/early/2017/12/01/120378.
%
% USAGE:
% >>[X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
% >>[X,genelist]=sc_selectg(X,genelist);
% >>figure; s1=run_phate(X,3,true,true,genelist);  % view genes
% >>figure; s2=run_phate(X,3,true);          % view cells
% >>figure; scatter3(s2(:,1),s2(:,2),s2(:,3),10,'filled'); % using output S

if nargin < 2, ndim = 3; end
if nargin < 3, plotit = false; end
if nargin < 4, bygene = false; end
if nargin < 5, genelist = []; end

pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    % phate and its helpers (svdpca, randPCA, knee_pt, ...) come off the
    % path again when this function returns.
    phatecleanup = pkg.i_addpathtemp( ...
        fullfile(fileparts(pw1), 'external', 'ml_PHATE'));   %#ok<NASGU>
end
% gene_names=cellstr(gl123);
% PHATE on data (rows: samples, columns: features)
% data=X';
% % library size normalization
% libsize = sum(data,2);
% data = bsxfun(@rdivide, data, libsize) * median(libsize);

if bygene
    data = sc_norm(X', 'type', 'libsize');
else
    data = sc_norm(X, 'type', 'libsize');
end

% The following transpose is necessary to make the input dim right.
data = data';

% sqrt transform
% data = log(data+1);
data = sqrt(data);


% The bundled external/ml_PHATE/randPCA.m does a bare "warning off" at
% line 103 with nothing to undo it, so one PHATE run left warnings
% disabled for the rest of the MATLAB session -- silencing every later
% warning the user relies on. The same defect in external/ml_SinNLRR was
% fixed this way in 0827873 and in external/ml_MAGIC alongside this
% change; +pkg/e_randPCA.m, the in-toolbox copy of that very file, has the
% offending line commented out already. Restoring here rather than editing
% the third-party file keeps the fix in code we own.
warnState = warning();
restoreWarn = onCleanup(@() warning(warnState));

% That same randPCA also opens its randomized branch with rng('default')
% and never puts the stream back, so one PHATE run left the whole session
% parked on seed 0 -- every later stochastic step, t-SNE, UMAP, k-means, a
% bootstrap or the user's own code, drawing from that instead of the seed
% they set. Restored here for the same reason as the warning state above:
% the defect is in third-party code, the fix belongs in code we own.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));

s = phate(data, 't', 20, 'ndim', ndim, 'k', 10, 'npca', min([100, size(X,2)]));

%%
if plotit
    if bygene
        [lgu, dropr, lgcv] = sc_genestat(X, [], false);
        colorby = 'mean';
        switch colorby
            case 'mean'
                C = lgu;
            case 'dropout'
                C = dropr;
            case 'cv'
                C = lgcv;
        end
    else
        C = sum(X, 1); % by library size
    end
    switch ndim
        case 2
            scatter(s(:, 1), s(:, 2), 10, C, 'filled');
            % colormap(jet)
            % set(gca,'xticklabel',[]);
            % set(gca,'yticklabel',[]);
            % axis tight
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            title 'PHATE'
        case 3
            scatter3(s(:, 1), s(:, 2), s(:, 3), 10, C, 'filled');
            xlabel 'PHATE1'
            ylabel 'PHATE2'
            zlabel 'PHATE3'
            title 'PHATE 3D'
    end
    if ~isempty(genelist)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1, genelist};
    end
end
end

function txt = i_myupdatefcn1(~, event_obj, g)
% Data-tip text: the gene under the cursor. Local because the shared
% copies live in private/ folders that +run cannot see.
txt = {g(event_obj.DataIndex)};
end
