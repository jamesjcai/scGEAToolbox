function callback_GlycoEnrich(src, ~)
%GUI.CALLBACK_GLYCOENRICH  Menu callback: enrichment/depletion of
%glyco-module expressing cells across cell groups.
%
%   Wraps GLY.ENRICH, which asks, for every glyco module and every group of
%   cells, whether that group holds more or fewer cells EXPRESSING the
%   module than the rest of the data does. The answer is a log2 odds ratio
%   from a 2-by-2 table, with a Fisher exact test and a BH-adjusted p-value.
%
%   The grouping is chosen rather than assumed: glycosylation programs are
%   remodelled between conditions as well as between cell types, and which
%   of the two is being asked about is the analysis, not a detail.
%
% See also GLY.ENRICH, GLY.DETECT, GUI.I_GETCELLGROUPS.

[FigureHandle, sce] = gui.gui_getfigsce(src);

[grp, groupBy] = gui.i_getcellgroups(sce, FigureHandle, Ask=true, ...
    Prompt="Group cells by:");
if isempty(grp), return; end

minDetected = in_askmindetected(FigureHandle);
if isempty(minDetected), return; end

fw = gui.myWaitbar(FigureHandle);
try
    [T, info] = gly.enrich(sce.X, sce.g, grp, MinDetected=minDetected);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.enrich');
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(T)
    gui.myWarndlg(FigureHandle, ['No glyco module could be tested on this ', ...
        'dataset. Too few of the collection''s genes are present.']);
    return;
end

in_warnaboutdepth(info, groupBy, FigureHandle);
gui.i_viewtable(T, FigureHandle);
gui.i_exporttable(T, true, 'Tglycoenrich', 'GlycoEnrichTable', ...
    [], [], FigureHandle);
end


function minDetected = in_askmindetected(parentfig)
% 1 reproduces the published definition (mean pathway expression nonzero).
% On a large module it makes almost every cell an expressing cell and the
% odds ratio goes flat, which is why it is exposed rather than fixed.

minDetected = [];
answer = gui.i_inputdlg(['Genes of a module that must be detected in a cell ', ...
    'for it to count as expressing (1 = the published definition):'], ...
    '1', parentfig);
if isempty(answer), return; end

v = str2double(answer{1});
if ~isfinite(v) || v < 1 || v ~= fix(v)
    gui.myErrordlg(parentfig, ...
        'Enter a whole number of 1 or more.', 'gly.enrich');
    return;
end
minDetected = v;
end


function in_warnaboutdepth(info, groupBy, parentfig)
% Expressing fraction is a detection statistic, so a group sequenced more
% deeply expresses more of everything. INFO.rhoDepth measures exactly that
% per module, and the whole table is worth less where it runs high - so say
% so before the user reads the table, not after.

if ~isfield(info, 'rhoDepth') || all(isnan(info.rhoDepth))
    if numel(info.groups) < 3
        gui.myHelpdlg(parentfig, [ ...
            "Only two groups, so the depth-confounding check (the correlation " + ...
            "between each module's expressing fraction and group sequencing " + ...
            "depth) could not be computed - it needs three or more."
            ""
            "Run the glycogene detection depth check on the same grouping to " + ...
            "see whether the groups are comparably sequenced."], ...
            'Glyco-module enrichment');
    end
    return;
end

confounded = sum(abs(info.rhoDepth) >= 0.8);
if confounded == 0, return; end

gui.myWarndlg(parentfig, sprintf( ...
    ['In %d of %s, the expressing fraction is correlated with ', ...
     'group sequencing depth at |rho| >= 0.8 (grouped by %s).\n\n', ...
     'That is the signature of a depth artefact rather than glycobiology: ', ...
     'a more deeply sequenced group expresses more of everything. Treat ', ...
     'those rows as unresolved, and check the depth balance with the ', ...
     'glycogene detection depth check before reporting them.'], ...
    confounded, pkg.i_plural(numel(info.rhoDepth), 'module'), groupBy));
end
