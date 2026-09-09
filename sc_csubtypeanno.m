function [sce] = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag)
%SC_CSUBTYPEANNO Annotate the subtypes of one primary cell type.
%
%   sce = sc_csubtypeanno(sce, cell_type_target)
%   sce = sc_csubtypeanno(sce, cell_type_target, formatid, speciestag)
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
%   formatid   0 subtype alone (default), 1 'Type_{Subtype}', 2 'Type (Subtype)'
%   speciestag 'human' (default) or 'mouse'
%
%   See also GUI.CALLBACK_SUBTYPEANNOTATION, PKG.I_MATCHPRIMARYTYPE,
%   PKG.E_DETERMINECELLTYPE.

if nargin < 4 || isempty(speciestag)
    speciestag = 'human';
end
speciestag = validatestring(lower(char(speciestag)), {'human', 'mouse'}, ...
    'sc_csubtypeanno', 'speciestag', 4);

if nargin < 3 || isempty(formatid)
    formatid = 0;
end

pw1 = fileparts(mfilename('fullpath'));
pth2 = fullfile(pw1, 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');

T = readtable(pth2);

if ~ismember(upper(cell_type_target), upper(string(unique(T.CellType))))
    error('Target cell type is not a primary cell type in cellsubtype.xlsx.');
end

% Every label that means the target type, not only the ones spelled the way the
% marker table spells it.
selectedidx = pkg.i_matchprimarytype(sce.c_cell_type_tx, cell_type_target) ...
    == string(cell_type_target);
selectedidx = reshape(selectedidx, size(sce.c_cell_type_tx));

if ~any(selectedidx)
    error('SCE.C_CELLTYPE_TXT does not contain the target cell type.');
end

idx = upper(string(T.CellType)) == upper(cell_type_target);

if ~any(idx)
    error('SC_CSUBTYPEANNO: Runtime Error.');
end

T = T(idx, :);

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
[pmarkerstr] = pkg.e_primarymarkers(cell_type_target, speciestag);
Tm = T(:, 2:3);
[Tm] = in_addprimarymarkers(Tm, pmarkerstr);

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

sce2 = copy(sce);
sce2 = sce2.selectcells(selectedidx); % OK

% Isolating one cell type leaves genes with no counts at all in it. They carry
% nothing for the re-embedding, and SC_SPLINEFIT drops them anyway - noisily,
% one warning per run. Dropping them here also takes them out of the marker
% scoring below, where a marker that is measured but never detected in this
% population was counted among the matched genes and diluted that subtype's
% score without adding to it.
sce2 = sce2.selectgenes(1, 1);   % keep genes with >=1 count in >=1 cell

sce2 = sce2.embedcells('tsne3d', true, true, 3);
sce2 = sce2.clustercells([], [], true);

[c, cL] = findgroups(string(sce2.c_cluster_id));

for ik = 1:max(c)
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
            a(k) = sprintf('%s_{%s}', b, a(k));
        end
    case 2
        for k = 1:length(a)
            a(k) = sprintf('%s (%s)', b, a(k));
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

