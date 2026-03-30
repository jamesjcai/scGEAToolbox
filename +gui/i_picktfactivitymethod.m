function [answer, id] = i_picktfactivitymethod(parentfig)
% I_PICKTFACTIVITYMETHOD  Select TF activity scoring method for sc_tfactivity.
%
% All methods use the DoRothEA TF→target database and differ only in how
% they aggregate target gene expression into a per-cell activity score.
% Listed in order of increasing computational cost (fastest → slowest).
%
% Methods 1–5 use positive regulatory relationships only (mor > 0).
% Method 6 (VIPER) uses the full signed regulon (activators + repressors).
%
% Display order  →  methodid  (ordered by computational complexity)
%   1. WMEAN (weighted mean) — matrix multiply           (fastest) [DEFAULT]
%   2. UCell                 — rank-based AUC
%   3. VIPER / aREA          — signed rank enrichment
%   4. ULM / decoupleR       — regression + p-values
%   5. Pseudo-inverse        — SVD-based least squares
%   6. NMF                   — iterative factorization  (slowest 🐢)

if nargin < 1, parentfig = []; end

items = { ...
    '① WMEAN                              — weighted mean of target genes (fastest)', ...
    '② UCell [PMID:34285779]             — rank-based AUC', ...
    '③ VIPER / aREA [PMID:27322546]      — signed regulon rank enrichment', ...
    '④ ULM / decoupleR                   — regression with p-values', ...
    '⑤ Pseudo-inverse                    — SVD-based least squares', ...
    '⑥ NMF [PMID:33135076]              — iterative factorization 🐢 (slowest)'};

% Default selection: AddModuleScore (index 1 in this list)
default_idx = 1;

if gui.i_isuifig(parentfig)
    [idx, ok] = gui.myListdlg(parentfig, items, ...
        'Select TF activity method (ordered fastest → slowest):', ...
        items(default_idx), true,  true, [440, 450]);
else
    [idx, ok] = listdlg( ...
        'PromptString', 'Select TF activity method (fastest → slowest):', ...
        'SelectionMode', 'single', ...
        'ListString', items, ...
        'ListSize', [440, 180], ...
        'InitialValue', default_idx);
end

if ~ok || isempty(idx)
    answer = [];
    id     = [];
    return;
end

answer = items{idx};
id = idx; % display position matches methodid in sc_tfactivity (1:1)
end
