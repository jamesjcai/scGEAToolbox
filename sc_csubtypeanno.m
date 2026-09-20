function [sce] = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag, fw, opts)
%SC_CSUBTYPEANNO Annotate the subtypes of one primary cell type.
%
%   sce = sc_csubtypeanno(sce, cell_type_target)
%   sce = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag)
%   sce = sc_csubtypeanno(..., fw, opts)
%
%   Cells annotated as cell_type_target are isolated, re-embedded and
%   re-clustered on their own, each of the new clusters is scored against the
%   subtype markers of that type in assets/PanglaoDB/cellsubtypes.xlsx, and the
%   subtype labels are merged back into sce.c_cell_type_tx in place. Cells of
%   every other type keep the labels they had.
%
%   cell_type_target is a primary type of the marker table, e.g. "T cells".
%   Which cells belong to it is decided by PKG.I_MATCHPRIMARYTYPE rather than
%   by an exact string comparison, so a gross annotation that spells the type
%   as "T cells_{3}", "CD8+ T cells" or "T memory cells" is picked up as well.
%
%   Cells whose label already names a subtype of the target - "T regulatory
%   cells", "Plasma cells", "Interneurons", the 37 names PKG.I_SUBTYPEOVERLAP
%   knows celltypes.xlsx and cellsubtypes.xlsx both cover - are recognized as
%   the target's cells but left out of the run and keep their labels. Pass
%   OPTS.CELLSELECTION to subdivide them anyway.
%
%   formatid   0 subtype alone (default), 1 'Type_{Subtype}', 2 'Type (Subtype)'
%   speciestag 'human' (default) or 'mouse'
%   fw         waitbar context, [] for none
%   opts       struct of overrides, any field may be omitted:
%
%     SubtypeMarkers    a two-column table - name, comma separated markers -
%                       to score against instead of the bundled subtype table.
%                       Supplying it lifts the restriction to the primary
%                       types cellsubtypes.xlsx covers: any cell type can be
%                       subdivided, into subtypes of the caller's choosing.
%     AddPrimaryMarkers append the primary markers of CELL_TYPE_TARGET to every
%                       subtype's list. Default true for the bundled table,
%                       where the subtype lists are short and assume them, and
%                       false for a supplied one, where they would dilute a
%                       list the caller already considers complete.
%     CellSelection     logical index of the cells to subdivide, for when the
%                       target is a label rather than a primary type.
%                       Default: whatever PKG.I_MATCHPRIMARYTYPE matches.
%
%   See also GUI.CALLBACK_SUBTYPEANNOTATION, PKG.I_MATCHPRIMARYTYPE,
%   PKG.I_SUBTYPEOVERLAP, PKG.E_DETERMINECELLTYPE.

if nargin < 6, opts = struct(); end
if nargin < 5, fw = []; end

% Presence of the field, not its contents, is what selects the custom path:
% an empty marker table is a caller's mistake and has to say so, not fall back
% to the bundled one and annotate against markers nobody asked for.
usecustom = isstruct(opts) && isscalar(opts) && isfield(opts, 'SubtypeMarkers');
customtable = [];
if usecustom, customtable = opts.SubtypeMarkers; end
addprimary = in_optfield(opts, 'AddPrimaryMarkers', ~usecustom);
if nargin < 4 || isempty(speciestag)
    speciestag = 'human';
end
speciestag = validatestring(lower(char(speciestag)), {'human', 'mouse'}, ...
    'sc_csubtypeanno', 'speciestag', 4);

if nargin < 3 || isempty(formatid)
    formatid = 0;
end

if usecustom
    % The caller brought its own subtypes, so there is no primary type to
    % check the target against and nothing to look up in the bundled table.
    Tm = in_customtable(customtable);
else
    pw1 = fileparts(mfilename('fullpath'));
    pth2 = fullfile(pw1, 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');

    T = readtable(pth2);

    if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType))))
        error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
    end

    idx = upper(string(T.CellType)) == upper(cell_type_target);

    if ~any(idx)
        error('SC_CSUBTYPEANNO: Runtime Error.');
    end

    T = T(idx, :);
    Tm = T(:, 2:3);
end

selectedidx = in_optfield(opts, 'CellSelection', []);
nresolved = 0;
if isempty(selectedidx)
    % Every label that means the target type, not only the ones spelled the way
    % the marker table spells it.
    selectedidx = pkg.i_matchprimarytype(sce.c_cell_type_tx, cell_type_target) ...
        == string(cell_type_target);

    % ...except the labels that already name a subtype of it. celltypes.xlsx
    % and cellsubtypes.xlsx overlap, so a primary annotation run can hand back
    % "T regulatory cells" or "Plasma cells" alongside "T cells" and "B cells".
    % Those cells belong to the target - PKG.I_SUBTYPEOVERLAP is how they were
    % recognized as its cells in the first place - but re-clustering them would
    % replace a subtype the markers already decided with whichever subtype
    % their new cluster happens to score highest for. Coarser, and no more
    % likely to be right. They keep the label they have.
    [oprimary, ~, isoverlap] = pkg.i_subtypeoverlap(sce.c_cell_type_tx);
    isresolved = reshape(isoverlap & oprimary == string(cell_type_target), ...
        size(selectedidx));
    nresolved = sum(selectedidx & isresolved);
    selectedidx = selectedidx & ~isresolved;
end
selectedidx = reshape(logical(selectedidx), size(sce.c_cell_type_tx));

if ~any(selectedidx)
    if nresolved > 0
        % Not the same complaint: the cells are there, they are just all
        % subtyped already, and the caller has nothing left to gain.
        error('sc_csubtypeanno:AlreadySubtyped', ['All %d %s in the data ' ...
            'already carry a subtype label. There is nothing left to ' ...
            'subdivide. Collapse the subtypes back to the cell type first ' ...
            'if you want them re-annotated.'], nresolved, cell_type_target);
    end
    error('SCE.C_CELLTYPE_TXT does not contain the target cell type.');
end

% SPECIESTAG is passed on. It was documented and defaulted but never read:
% the local loader this replaces took only (pw1, cell_type_target) and
% loaded marker_hs.mat unconditionally, so a mouse dataset had HUMAN
% primary markers appended to every subtype's list. The two files really
% do differ -- marker_mm.mat sits beside marker_hs.mat and is not a copy.
% Over the primary types the subtype table covers, the shared fraction of
% each marker set runs from 0.73 to 1.00: Macrophages has 128 human
% markers against 147 mouse with 123 in common, Dendritic cells 121
% against 108 with 97. There are also 4 human-only and 8 mouse-only
% primary types.
if addprimary
    [pmarkerstr] = pkg.e_primarymarkers(cell_type_target, speciestag);
    [Tm] = in_addprimarymarkers(Tm, pmarkerstr);
end

Tw = pkg.e_markerweight(Tm);

wvalu = Tw.Var2;
wgene = string(Tw.Var1);
celltypev = string(Tm.SubType);

% UPPER, to match what it is compared against. PKG.E_DETERMINECELLTYPE
% splits MARKERGENEV and tests each gene with == against WGENE, which
% PKG.E_MARKERWEIGHT has already upper-cased, and against upper(sce.g). A
% marker arriving in mixed case therefore matches neither and is silently
% dropped from that subtype's score.
%
% This changes nothing for the shipped table: all 956 distinct markers in
% assets/PanglaoDB/cellsubtypes.xlsx are already upper case, as are all
% 1529 in the primary marker files appended above. It is here so that an
% edited or extended marker table behaves, rather than losing rows without
% saying so.
markergenev = upper(string(Tm.PositiveMarkers));

if ~isempty(fw)
    gui.myWaitbar(fw.FigureHandle, fw.fw, false, '', ...
        'Extracting major cell type for subtype annotation...', 0.1);
end
sce2 = copy(sce);
sce2 = sce2.selectcells(selectedidx); % OK

% Isolating one cell type leaves genes with no counts at all in it. They carry
% nothing for the re-embedding, and SC_SPLINEFIT drops them anyway - noisily,
% one warning per run. Dropping them here also takes them out of the marker
% scoring below, where a marker that is measured but never detected in this
% population was counted among the matched genes and diluted that subtype's
% score without adding to it.
sce2 = sce2.selectgenes(1, 1);   % keep genes with >=1 count in >=1 cell

if ~isempty(fw)
    gui.myWaitbar(fw.FigureHandle, fw.fw, false, '', ...
        'Reembedding extracted cells...', 0.2);
end
sce2 = sce2.embedcells('tsne3d', true, true, 3);

if ~isempty(fw)
    gui.myWaitbar(fw.FigureHandle, fw.fw, false, '', ...
        'Clustering extract cells...', 0.3);
end
% Ask for a cluster count instead of taking the default. CLUSTERCELLS
% defaults to round(NumCells/100, -1), which rounds to the nearest TEN and so
% is 0 -- raised to 1 -- for every population under 500 cells. One cluster
% means one score and one label, so the subtypes of a 300-cell population all
% came back as whichever subtype scored highest overall. An isolated cell type
% is a few hundred cells often enough that this was the usual outcome.
%
% At least one cluster per subtype, so each has somewhere to land, and up to
% two, so a subtype that splits is not forced back together; capped at one
% cluster per 25 cells, below which a cluster is too small to score.
kclust = in_clustercount(sce2.NumCells, height(Tm));
sce2 = sce2.clustercells(kclust, [], true);

[c, cL] = findgroups(string(sce2.c_cluster_id));

for ik = 1:max(c)
    if ~isempty(fw)
        gui.myWaitbar(fw.FigureHandle, fw.fw, false, '', ...
            sprintf('Annotating cell types using PanglaoDB for cluster %d of %d...', ik, max(c)),...
            0.3 + 0.7 * (ik/max(c)));
    end
    ptsSelected = c == ik;
    [Tct] = pkg.e_determinecelltype(sce2, ptsSelected, wvalu, ...
            wgene, celltypev, markergenev);

    ctxt = Tct.C1_Cell_Type{1};
    cL{ik} = ctxt;
end

sce2.c_cell_type_tx = string(cL(c));

sce.c_cell_type_tx(selectedidx) = in_formatsubtype(sce2.c_cell_type_tx, ...
    cell_type_target, formatid);
end

function [a] = in_formatsubtype(a, b, formatid)
switch formatid
    case 0
        return;
    case 1
        for k = 1:length(a)
            a(k) = sprintf('%s (%s)', b, a(k));
        end
    case 2
        for k = 1:length(a)
            a(k) = sprintf('%s_{%s}', b, a(k));
        end
    otherwise
        error('sc_csubtypeanno:InvalidFormat', 'Unknown formatid %d. Use 0, 1, or 2.', formatid);
end
end


function [Tm] = in_addprimarymarkers(Tm, pmarkerstr)
for k = 1:size(Tm, 1)
    a = string(Tm.PositiveMarkers{k});
    a = strtrim(a);
    a = erase(a, " ");
    a = strip(a, 'right', ',');
    Tm.PositiveMarkers{k} = char(string(a) + "," + pmarkerstr);
end
end


function [v] = in_optfield(opts, name, default)
v = default;
if isstruct(opts) && isscalar(opts) && isfield(opts, name) && ...
        ~isempty(opts.(name))
    v = opts.(name);
end
end

function [Tm] = in_customtable(T)
% A caller-supplied marker table, normalized to the two columns the rest of
% this function reads: SubType and PositiveMarkers. Accepts either name for
% either column, or a plain two-column table as PKG.I_PARSEMARKERLIST returns
% one.

if ~istable(T) || width(T) < 2 || isempty(T)
    error('sc_csubtypeanno:BadMarkerTable', ['SubtypeMarkers must be a ' ...
        'non-empty table with a name column and a marker gene column.']);
end

vars = string(T.Properties.VariableNames);
namecol = find(ismember(lower(vars), ["subtype", "celltype", "name", "var1"]), 1);
genecol = find(ismember(lower(vars), ...
    ["positivemarkers", "markers", "genes", "var2"]), 1);
if isempty(namecol), namecol = 1; end
if isempty(genecol) || genecol == namecol
    genecol = find((1:width(T)) ~= namecol, 1);
end

SubType = cellstr(strtrim(string(T{:, namecol})));
PositiveMarkers = cellstr(upper(erase(strtrim(string(T{:, genecol})), " ")));
Tm = table(SubType, PositiveMarkers);
end

function [k] = in_clustercount(ncells, nsubtypes)
k = max(2, min(2*nsubtypes, floor(ncells/25)));
k = max(1, min(k, ncells));
end
