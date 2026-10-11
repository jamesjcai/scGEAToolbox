function sce = qcfilterwhitelist(sce, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff, whitelist)
if nargin < 6, whitelist = []; end
if nargin < 5 || isempty(gnnumcutoff), gnnumcutoff = 200; end
if nargin < 4 || isempty(min_cells_nonzero), min_cells_nonzero = 15; end
if nargin < 3 || isempty(mtratio), mtratio = 0.15; end
if nargin < 2 || isempty(libszcutoff), libszcutoff = 1000; end

before = [sce.NumCells, sce.NumGenes];
if ~isempty(whitelist)
    assert(all(ismember(whitelist, sce.g)), ...
        'All whitelist genes must be present in sce.g.');
    [~, idxx] = ismember(whitelist, sce.g);
    Xresv = sce.X(idxx, :);
    % Capture the whitelist rows' gene attributes now, so they can be put
    % back alongside Xresv after the whitelist genes are re-appended below.
    attrResv = i_getGeneAttributeRows(sce, idxx);
end

[~, keptg, keptidxv] = sc_qcfilter(sce.X, sce.g, libszcutoff, mtratio, ...
min_cells_nonzero, gnnumcutoff);

for k = 1:length(keptidxv)
    sce = selectcells(sce, keptidxv{k}); % OK
end

[keptmask] = ismember(sce.g, keptg);
sce = sce.selectgenesbyindex(keptmask);
% Default methodid renames duplicates rather than dropping them, so the
% gene count - and therefore the gene attributes - are unchanged here.
[sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g);

if ~isempty(whitelist)
    for k = 1:length(keptidxv)
        Xresv = Xresv(:, keptidxv{k});
    end
    [found, idxx] = ismember(whitelist, sce.g);

    if any(found)
        % IDXX is subscripts from ISMEMBER, so the complement is built
        % rather than negated. SELECTGENESBYINDEX carries the gene
        % attributes across, which used to need the separate call below.
        idxx = idxx(found);
        keepmask = true(sce.NumGenes, 1);
        keepmask(idxx) = false;
        sce = sce.selectgenesbyindex(keepmask);
    end
    sce.X = [sce.X; Xresv];
    sce.g = [sce.g; whitelist(:)];
    sce = i_appendGeneAttributeRows(sce, attrResv);
end

% Stored, not printed: this also runs inside other analyses. Entry points
% print it (PKG.I_METHODSSHOW).
pkg.i_methodsreport(sce, "qc", pkg.i_qcmethodstext(LibSize=libszcutoff, ...
    MtRatio=mtratio, MinCells=min_cells_nonzero, MinGenes=gnnumcutoff, ...
    NumWhitelist=numel(whitelist), CellsBefore=before(1), ...
    GenesBefore=before(2), CellsAfter=sce.NumCells, ...
    GenesAfter=sce.NumGenes), "cai2020", Print=false);
end
