function sce = qcfilterwhitelist(sce, libszcutoff, mtratio, ...
    min_cells_nonzero, gnnumcutoff, whitelist)
if nargin < 6, whitelist = []; end
if nargin < 5 || isempty(gnnumcutoff), gnnumcutoff = 200; end
if nargin < 4 || isempty(min_cells_nonzero), min_cells_nonzero = 15; end
if nargin < 3 || isempty(mtratio), mtratio = 0.15; end
if nargin < 2 || isempty(libszcutoff), libszcutoff = 1000; end

% if issparse(obj.X)
%     try
%         obj.X = full(obj.X);
%     catch
%     end
% end

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
sce.X = sce.X(keptmask, :);
sce.g = sce.g(keptmask);
sce = i_filterGeneAttributes(sce, keptmask);
% Default methodid renames duplicates rather than dropping them, so the
% gene count - and therefore the gene attributes - are unchanged here.
[sce.X, sce.g] = sc_rmdugenes(sce.X, sce.g);

if ~isempty(whitelist)
    for k = 1:length(keptidxv)
        Xresv = Xresv(:, keptidxv{k});
    end
    [found, idxx] = ismember(whitelist, sce.g);

    if any(found)
        idxx = idxx(found);
        keepmask = true(sce.NumGenes, 1);
        keepmask(idxx) = false;
        sce.X(idxx, :) = [];
        sce.g(idxx) = [];
        sce = i_filterGeneAttributes(sce, keepmask);
    end
    sce.X = [sce.X; Xresv];
    sce.g = [sce.g; whitelist(:)];
    sce = i_appendGeneAttributeRows(sce, attrResv);
end

end
