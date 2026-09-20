function callback_GlycoCcc(src, ~)
%GUI.CALLBACK_GLYCOCCC  Menu callback: glyco-lectin cell-cell communication.
%
%   Wraps GLY.CCC. Unlike a protein ligand-receptor engine, the "ligand"
%   here is a glycan determinant, which RNA-seq cannot see: it is stood in
%   for by the sender's biosynthetic capacity (a glyco-module score from
%   GLY.STATE), while the "receptor" is a lectin gene actually measured on
%   the receiver. Each sender to receiver edge is scored as mean sender
%   module score times mean receiver lectin expression, against a
%   permutation null.
%
%   This is a different channel from GUI.CALLBACK_GLYCOWEIGHT, not an
%   earlier stage of it: this one finds glycan-mediated edges, that one
%   discounts protein ligand-receptor edges by their glycan context. The
%   two tables have different columns and neither feeds the other.
%
% See also GLY.CCC, GLY.LECTINMAP, GLY.STATE, GUI.CALLBACK_GLYCOWEIGHT.

[FigureHandle, sce] = gui.gui_getfigsce(src);

[labels, ~] = gui.i_getcellgroups(sce, FigureHandle);
if isempty(labels), return; end

% GLY.CCC multiplies a module score by a mean expression level, so it wants
% normalized input. On raw counts the lectin side is dominated by library
% size and the edge ranking follows depth.
if in_looksrawcounts(sce.X)
    answer = gui.myQuestdlg(FigureHandle, [ ...
        "The counts in this dataset look like raw integers."
        ""
        "This analysis multiplies a glyco-module score by mean lectin " + ...
        "expression, which wants normalized (e.g. log1p) input. On raw " + ...
        "counts the receiver side follows library size, so the edge " + ...
        "ranking will partly be a depth ranking."
        ""
        "Continue anyway?"]);
    if ~strcmp(answer, 'Yes'), return; end
end

[~, methodid] = gui.i_pickscoremethod([], FigureHandle);
if isempty(methodid), return; end

nPerm = in_asknumber(FigureHandle, ['Permutations for the null ', ...
    '(more is slower and resolves smaller p-values):'], '100', 10, 'gly.ccc');
if isempty(nPerm), return; end

fw = gui.myWaitbar(FigureHandle);
try
    result = gly.ccc(sce.X, sce.g, labels, methodid=methodid, n_perm=nPerm);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'gly.ccc');
    return;
end
gui.myWaitbar(FigureHandle, fw);

T = result.T_interactions;
if isempty(T)
    gui.myHelpdlg(FigureHandle, [ ...
        "No glycan-lectin edge reached significance on this dataset."
        ""
        "With " + string(nPerm) + " permutations the smallest reachable " + ...
        "p-value is " + string(1/nPerm) + ", so a null result at a low " + ...
        "permutation count may only mean the test could not resolve one. " + ...
        "Check that the cognate lectin genes are detected here before " + ...
        "reading this as an absence of glycan-mediated signalling."], ...
        'Glyco-lectin communication');
    return;
end

gui.myHelpdlg(FigureHandle, [ ...
    sprintf('%d significant edge(s) across %d cell types, from %d epitope row(s).', ...
        height(T), numel(result.cell_types), numel(result.epitopes))
    ""
    "The sender side is biosynthetic capacity, not a measured glycan. An " + ...
    "edge says the sender can build the determinant and the receiver " + ...
    "expresses a lectin that reads it - not that the glycan was observed."], ...
    'Glyco-lectin communication');

gui.i_viewtable(T, FigureHandle);
gui.i_exporttable(T, true, 'Tglycoccc', 'GlycoCccTable', [], [], FigureHandle);
end


function tf = in_looksrawcounts(X)
% Integrality on a sample of the nonzeros. Normalized or log1p data is
% almost never integral, and sampling keeps this cheap on a large matrix.

v = nonzeros(X);
if isempty(v), tf = false; return; end
if numel(v) > 2000
    v = v(round(linspace(1, numel(v), 2000)));
end
v = full(double(v));
tf = all(v == fix(v));
end


function v = in_asknumber(parentfig, prompt, definput, minimum, errid)
% One whole number, validated. Returns empty on cancel or a bad entry.

v = [];
answer = gui.i_inputdlg(prompt, definput, parentfig);
if isempty(answer), return; end
x = str2double(answer{1});
if ~isfinite(x) || x < minimum || x ~= fix(x)
    gui.myErrordlg(parentfig, sprintf( ...
        'Enter a whole number of %d or more.', minimum), errid);
    return;
end
v = x;
end
