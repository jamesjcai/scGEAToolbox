function [sce] = sc_mergesces(sces, method, keepbatchid, forcenoappendix)
% Merges two SCE objects
% Usage: [sce]=sc_mergesces({sce1,sce2},'intersect');
% See also: SC_MERGEDATA

if nargin < 4, forcenoappendix = false; end
if nargin < 3
    keepbatchid = true;
end

if nargin < 2 || isempty(method)
    method = "intersect";
end

method = string(method);                 % normalize type
validMethods = ["intersect","union"];
method = validatestring(method, validMethods);
method = string(method);                 % convert back to string

if ~iscell(sces), error('SCES is not a cell array.'); end
if length(sces) < 2, error('At least two SCE required.'); end
for k = 1:length(sces)
    if ~isa(sces{k}, 'SingleCellExperiment')
        error('sces{%d} is not a SingleCellExperiment object.', k);
    end
end

needappendix = false;
sce = sces{1};
ccell = cell(length(sces), 1);
ccell{1} = ones(sce.NumCells, 1);
for k = 2:length(sces)
    ccell{k} = k * ones(sces{k}.NumCells, 1);
    [sce, hasidoverlapx] = i_merge2sces(sce, sces{k}, method);
    if hasidoverlapx
        needappendix = true;
    end
end
c = vertcat(ccell{:});
if ~keepbatchid || isscalar(unique(sce.c_batch_id))
    sce.c_batch_id = c;
end
if ~forcenoappendix
    if needappendix
        sce.c_batch_id = strcat(string(sce.c_batch_id), "_", string(c));
        disp('A suffix is added to SCE.C_BATCH_ID to distinguish cells'' original batch IDs.');
    end
end

end


function [sce, hasidoverlap] = i_merge2sces(sce1, sce2, method)
hasidoverlap = false;
if nargin < 3, method = 'intersect'; end
[X, g, ~] = sc_mergedata(sce1.X, sce2.X, ...
    sce1.g, sce2.g, method);
sce = SingleCellExperiment(X, g);
sce.metadata = [sce1.metadata; sce2.metadata];
sce.c = [sce1.c; sce2.c];
try
    sce.s = [sce1.s; sce2.s];
catch ME
    warning('sc_mergesces:EmbeddingMerge', ...
        'Could not merge embeddings (%s). Using random s.', ME.message);
    sce.s = randn(size(X, 2), 3);
end
% sce.c_batch_id=c;

if ~isempty(sce1.c_batch_id) && ~isempty(sce2.c_batch_id)
    if ~isstring(sce1.c_batch_id)
        sce1.c_batch_id = string(sce1.c_batch_id);
    end
    if ~isstring(sce2.c_batch_id)
        sce2.c_batch_id = string(sce2.c_batch_id);
    end
    sce.c_batch_id = [sce1.c_batch_id; sce2.c_batch_id];
    % intersect(sce1.c_batch_id, sce2.c_batch_id)
    % pause
    hasidoverlap = ~isempty(intersect(sce1.c_batch_id, sce2.c_batch_id));
end

if ~isempty(sce1.c_cell_cycle_tx) && ~isempty(sce2.c_cell_cycle_tx)
    sce.c_cell_cycle_tx = [sce1.c_cell_cycle_tx; sce2.c_cell_cycle_tx];
end

if ~isempty(sce1.c_cell_type_tx) && ~isempty(sce2.c_cell_type_tx)
    if ~isstring(sce1.c_cell_type_tx)
        sce1.c_cell_type_tx = string(sce1.c_cell_type_tx);
    end
    if ~isstring(sce2.c_cell_type_tx)
        sce2.c_cell_type_tx = string(sce2.c_cell_type_tx);
    end
    sce.c_cell_type_tx = [i_remove_affix(sce1.c_cell_type_tx); ...
        i_remove_affix(sce2.c_cell_type_tx)];
end

if ~isempty(sce1.c_cluster_id) && ~isempty(sce2.c_cluster_id)
    sce.c_cluster_id = [sce1.c_cluster_id; sce2.c_cluster_id];
end

if ~isempty(sce1.c_cell_id) && ~isempty(sce2.c_cell_id)
    sce.c_cell_id = [sce1.c_cell_id; sce2.c_cell_id];
end

try
    if ~isempty(sce1.list_cell_attributes) && ~isempty(sce2.list_cell_attributes)
        for k = 1:2:min([length(sce1.list_cell_attributes), length(sce2.list_cell_attributes)])
            if strcmp(sce1.list_cell_attributes{k}, ...
                    sce2.list_cell_attributes{k})
                sce.list_cell_attributes{k} = sce1.list_cell_attributes{k};
                sce.list_cell_attributes{k+1} = [sce1.list_cell_attributes{k+1}; ...
                    sce2.list_cell_attributes{k+1}];
            end
        end
    end
catch ME
    warning('sc_mergesces:CellAttrMerge', ...
        'SCE.LIST_CELL_ATTRIBUTES not merged: %s', ME.message);
end

try
    a = fieldnames(sce.struct_cell_clusterings);
    for k = 1:length(a)
        c1 = sce1.struct_cell_clusterings.(a{k});
        c2 = sce2.struct_cell_clusterings.(a{k});
        if ~isempty(c1) && ~isempty(c2)
            sce.struct_cell_clusterings.(a{k}) = [c1; c2];
        end
    end
catch ME
    warning('sc_mergesces:ClusteringMerge', ...
        'SCE.STRUCT_CELL_CLUSTERINGS not merged: %s', ME.message);
end

try
    a = fieldnames(sce.struct_cell_embeddings);
    for k = 1:length(a)
        c1 = sce1.struct_cell_embeddings.(a{k});
        c2 = sce2.struct_cell_embeddings.(a{k});
        if ~isempty(c1) && ~isempty(c2)
            sce.struct_cell_embeddings.(a{k}) = [c1; c2];
        end
    end
catch ME
    warning('sc_mergesces:EmbeddingStructMerge', ...
        'SCE.STRUCT_CELL_EMBEDDINGS not merged: %s', ME.message);
end

% Modalities are merged rather than dropped. Unlike the blocks above this one
% is NOT wrapped in try/catch: a modality that silently disappears on merge is
% exactly the failure the PRESENT mask exists to prevent, and swallowing the
% error would hide it.
sce = i_mergemodalities(sce, sce1, sce2, method);

end


function sce = i_mergemodalities(sce, sce1, sce2, method)
%I_MERGEMODALITIES Carry both inputs' modalities onto the merged cell axis.
%
%   Three cases, one per modality name in the union of the two inputs:
%
%   In BOTH  - the feature axes are merged with the SAME method as the genes,
%              via SC_MERGEDATA (which is generic over row names), and the
%              cells concatenate. An 'intersect' merge therefore keeps the
%              peaks or antibodies common to both panels, exactly as it keeps
%              the genes common to both.
%
%   In ONE   - the modality still belongs on the merged object, because its
%              cells are half of the merged cells. The other half gets zero
%              columns and a FALSE entry in PRESENT. That is the same
%              distinction the loader draws for a cell a modality never
%              measured: absent is not a measured zero, and collapsing the two
%              would claim an assay covered cells it never saw.
%
%   Feature annotations survive only when the merge leaves that modality's
%   feature axis unchanged; after an intersect or union they describe features
%   that may no longer be there, so they are dropped with a warning rather
%   than silently re-indexed onto the wrong rows.

names = union(sce1.modalityNames(), sce2.modalityNames(), 'stable');
if isempty(names), return, end

n1 = sce1.NumCells;
n2 = sce2.NumCells;

for k = 1:numel(names)
    name = names(k);
    m1 = sce1.getModality(name);
    m2 = sce2.getModality(name);

    if ~isempty(m1) && ~isempty(m2)
        [X, f] = sc_mergedata(m1.X, m2.X, m1.f, m2.f, method);
        present = [m1.present, m2.present];
        type = m1.type;
        if m1.type ~= m2.type
            warning('sc_mergesces:ModalityTypeMismatch', ...
                ['Modality "%s" is type "%s" in one input and "%s" in the ', ...
                 'other; keeping "%s".'], name, m1.type, m2.type, m1.type);
        end
        source = m1;
        if isempty(f)
            warning('sc_mergesces:ModalityNoCommonFeatures', ...
                ['Modality "%s" has no features in common between the two ', ...
                 'inputs under ''%s''; it is attached with zero features.'], ...
                name, method);
        end
    elseif ~isempty(m1)
        [X, f] = i_padmodality(m1, 0, n2);
        present = [m1.present, false(1, n2)];
        type = m1.type;
        source = m1;
        warning('sc_mergesces:ModalityPartialCoverage', ...
            ['Modality "%s" exists in only one input, so %d of %d merged ', ...
             'cells are marked absent in its PRESENT mask.'], ...
            name, n2, n1 + n2);
    else
        [X, f] = i_padmodality(m2, n1, 0);
        present = [false(1, n1), m2.present];
        type = m2.type;
        source = m2;
        warning('sc_mergesces:ModalityPartialCoverage', ...
            ['Modality "%s" exists in only one input, so %d of %d merged ', ...
             'cells are marked absent in its PRESENT mask.'], ...
            name, n1, n1 + n2);
    end

    sce.setModality(name, X, f, type, present);

    % Re-attach annotations only if this modality's features came through
    % untouched; otherwise row j no longer describes feature j.
    if width(source.featureAnn) > 0
        if isequal(string(f(:)), string(source.f(:)))
            for v = string(source.featureAnn.Properties.VariableNames)
                sce.setFeatureAnnotation(name, v, source.featureAnn.(v));
            end
        else
            warning('sc_mergesces:ModalityAnnotationDropped', ...
                ['Modality "%s" changed its feature axis on merge, so its ', ...
                 '%d feature annotation(s) were dropped.'], ...
                name, width(source.featureAnn));
        end
    end
end
end


function [X, f] = i_padmodality(m, nBefore, nAfter)
% Place a modality's cells inside a wider axis, zero-filling the rest.
nf = numel(m.f);
if issparse(m.X)
    X = [sparse(nf, nBefore), m.X, sparse(nf, nAfter)];
else
    X = [zeros(nf, nBefore, 'like', m.X), m.X, zeros(nf, nAfter, 'like', m.X)];
end
f = m.f;
end


function b = i_remove_affix(a)
b = a;
for k = 1:length(a)
    idx = strfind(a(k), '_{');
    if ~isempty(idx) && idx(1) > 0
        b(k) = extractBefore(a(k), idx(1));
    end
end
end
