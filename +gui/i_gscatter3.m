function [h] = i_gscatter3(s, c, methodid, targetc, hAx, opts)
%I_GSCATTER3  Colour an embedding by a grouping variable or by a value.
%
%   h = gui.i_gscatter3(s, c)
%   h = gui.i_gscatter3(s, c, methodid, targetc, hAx)
%   h = gui.i_gscatter3(___, Continuous=true)
%
%   S is NumCells-by-2 or -by-3 embedding coordinates. C is one value per
%   cell. METHODID picks the drawing style: 1 SCATTER (default), 2
%   GSCATTER, 3 SCATTER with every cell outside TARGETC faded out. HAX is
%   the axes to draw into, [] for the current one.
%
%   C is treated one of two ways.
%
%   A GROUPING VARIABLE - cluster ids, batch ids, cell types, a logical
%   mask - is reduced by FINDGROUPS and drawn with one colour per level.
%   The colour index is the level's rank, which callers rely on: SCGEATOOL
%   labels colour index k with its own cL(k), correct only while this
%   function keeps that order. FINDGROUPS is applied to the NUMBERS, not to
%   STRING(c): routing numeric C through STRING() ordered the groups
%   lexicographically, so from ten groups on the colours were scrambled
%   against the values - with clusters 1..12, cluster 10 took colour index
%   2 and cluster 2 took index 5.
%
%   A CONTINUOUS VALUE - a likelihood, a pattern weight, a pseudotime, a
%   percentage - is handed to SCATTER as CData unchanged, so CLim is the
%   range of the data and a plain COLORBAR reports the values themselves.
%   Reducing one of these by FINDGROUPS still drew a usable picture, since
%   rank is monotone in value, but it set CLim to [1, number of distinct
%   values]: a MELD likelihood running 0.45 to 0.53 came out labelled 1 to
%   295, and a CoGAPS pattern weight of 0.00001 to 0.0148 likewise. The
%   numbers read as cell counts and could not be compared against a second
%   figure, or against a threshold.
%
%   CONTINUOUS decides which. The default, "auto", calls C continuous when
%   it is numeric and holds any non-integer value. That is deliberately
%   one-sided: every identifier in this toolbox is integral, so no existing
%   grouping can be reclassified by it, while an integer-valued measurement
%   - a library size, a detected-gene count, a raw UMI count - is left as a
%   grouping unless the caller says otherwise. Pass Continuous=true for
%   those, or Continuous=false to force levels.
%
%   Only METHODID 1 can carry values: GSCATTER draws one series per level
%   and method 3 selects the level equal to TARGETC, so both are defined in
%   terms of groups and C is reduced for them whatever CONTINUOUS says.
%
%   NaN becomes one more level on the grouping path, placed last, so those
%   cells stay on screen. On the continuous path they do not: SCATTER draws
%   no marker for a CData of NaN, and there is no colour on a continuous
%   axis that means "no value". Drawing them separately in grey was tried
%   and gives up more than it buys, because H would then hold fewer points
%   than there are cells - and H is indexed by cell elsewhere. SCGEATOOL
%   overwrites app.h.CData with a full-length vector, its brush maps the
%   selected point indices straight to cells, and GUI.I_CELLTYPEDATATIP
%   builds one datatip row per point. Row alignment is worth more than the
%   grey markers, so a caller that needs NaN cells shown has to drop them
%   from S itself.
%
%   H holds one point per row of S on both paths.
%
% See also GUI.I_GSCATTER3B, PKG.I_MYCOLORLINES, GUI.I_CELLTYPEDATATIP.

arguments
    s
    c = []
    methodid = 1
    targetc = 1
    hAx = []
    opts.Continuous = "auto"
end

if isempty(c), c = ones(size(s, 1), 1); end
if isempty(methodid), methodid = 1; end

x = s(:, 1);
y = s(:, 2);
if size(s, 2) >= 3, z = s(:, 3); end

iscont = in_iscontinuous(c, opts.Continuous) && methodid == 1;

% STRING() throws on a sparse array ("Conversion to string from sparse
% single is not possible"), and SCE.X is single sparse from R2025a on, so
% anything derived from it by summing - library size, detected genes, the
% mt ratio - arrives here sparse. C is also handed to SCATTER as CData
% below, which rejects sparse too. Normalise once here rather than
% trusting every caller: a colour variable is never usefully sparse.
% Non-numeric C (cellstr, string, categorical) is left alone for
% FINDGROUPS, and is never continuous.
if isnumeric(c) || islogical(c)
    c = full(double(c(:)));
end

if iscont
    h = in_scattervalues(x, y, s, c, hAx);
    return;
end

if isnumeric(c)
    c = findgroups(c);
    % FINDGROUPS returns NaN for NaN input. Under the old STRING() path
    % those became an ordinary "NaN" group with a colour of its own;
    % leaving them NaN would instead drop the points from the plot without
    % saying so, so keep them as one group, placed last.
    isn = isnan(c);
    if any(isn)
        c(isn) = max([0; c(~isn)]) + 1;
    end
else
    c = findgroups(string(c));
end
kc = numel(unique(c));

switch methodid
    case 1
        if size(s, 2) == 2
            if isempty(hAx)
                h = scatter(x, y, 10, c);
            else
                h = scatter(hAx, x, y, 10, c);
            end
        elseif size(s, 2) >= 3

            if isempty(hAx)
                h = scatter3(x, y, z, 10, c);
            else
                assert(numel(c)==numel(x));
                h = scatter3(hAx, x, y, z, 10, c);
            end
        end
    case 2
        if size(s, 2) == 2
            h = gscatter(x, y, c, [], [], 5, 'off');
        elseif size(s, 2) >= 3
            h = gscatter3b(x, y, z, c, [], [], 5, 'off');
        end
        box off
    case 3
        h = gui.i_gscatter3(s, c, 1);
        h.MarkerEdgeAlpha = 0;
        hold on
        idx = c == targetc;
        h = gui.i_gscatter3(s(idx, :), c(idx), 1);
        hold off
end


if isempty(hAx)
    grid on
    colormap(pkg.i_mycolorlines(kc));
else
    grid(hAx, 'on');
    colormap(hAx, pkg.i_mycolorlines(kc));
end

end


function tf = in_iscontinuous(c, setting)
% Whether C carries values rather than levels. The "auto" rule is one-sided
% on purpose - see the help above.

if ~(isnumeric(c) || islogical(c))
    tf = false;                      % text and categorical are always levels
    return;
end

if ~(isstring(setting) || ischar(setting))
    tf = logical(setting);
    return;
end
if ~strcmpi(setting, "auto")
    error("gui:i_gscatter3:badContinuous", ...
        'Continuous must be true, false, or "auto".');
end

v = full(double(c(:)));
v = v(isfinite(v));
tf = ~isempty(v) && any(v ~= round(v));
end


function h = in_scattervalues(x, y, s, c, hAx)
% The continuous path: CData as given, so CLim is the data's own range and
% a plain COLORBAR reports values. Every row of S is plotted, NaN included,
% because H is indexed by cell elsewhere - see the note in the help.

if isempty(hAx), hAx = gca(); end

if size(s, 2) >= 3
    h = scatter3(hAx, x, y, s(:, 3), 10, c);
else
    h = scatter(hAx, x, y, 10, c);
end

% 256 levels of the same map the grouping path reaches for above seven
% groups, so a figure does not change colour scheme when it changes scale.
grid(hAx, 'on');
colormap(hAx, turbo(256));
end
