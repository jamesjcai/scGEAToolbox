function [labels, report, percell, percellweight] = i_pooledannotate(gid, stage1, annotatecells, opts)
%I_POOLEDANNOTATE Settle a pooled annotation, labelling cells only where the pool is unclear.
%
%   labels = pkg.i_pooledannotate(gid, stage1, annotatecells)
%   [labels, report, percell, percellweight] = ...
%       pkg.i_pooledannotate(___, MinConfidence=0.5)
%
%   Stage one has already happened when this is called: the caller built one
%   pseudobulk profile per pool - and one per half of each pool - and sent
%   them to the classifier, so STAGE1 carries a label and a confidence for
%   every pool without a single cell having been classified on its own. This
%   decides which of those pools to believe, labels the cells of the rest,
%   and assembles the answer.
%
%   A pool is settled on its pseudobulk when its confidence clears
%   MinConfidence and its two halves were given the same label as the whole.
%   That pair of tests is deliberate, and neither is redundant:
%
%     - Confidence says the pooled profile landed somewhere the reference
%       atlas agrees about. It says nothing about whether the pool is one
%       cell type, because a pseudobulk of two types is a profile in its own
%       right and can land confidently on a third.
%     - Halves agreeing says the pool is homogeneous. Splitting along the
%       pool's own principal direction is what separates two populations if
%       there are two; if there is only one, the split is arbitrary and both
%       halves reproduce the whole. This is a local stand-in for
%       SCimilarity's query_coherence, which k-means a cluster into
%       sub-centroids and measures how far their reference neighbourhoods
%       overlap the whole cluster's - the same question, without needing the
%       21 GB cell-search index that query_coherence queries.
%
%   A pool failing either test is labelled cell by cell, in one batched
%   call, and then settled on the share its winning label holds. That vote
%   is weighted by whatever confidence the classifier returned with those
%   labels, so a handful of cells it was sure about outweighs more cells it
%   was not. A pool that still has no majority keeps its cells' own labels
%   and is reported as mixed, rather than being flattened to a consensus
%   untrue of most of it.
%
%   WHY A PSEUDOBULK AND NOT A SAMPLE. An earlier version of this function
%   labelled 5 cells per pool and took the majority. The pseudobulk is
%   cheaper - three profiles per 50-cell pool against five cells - and uses
%   every cell's counts instead of a tenth of them. SCimilarity builds it
%   the same way in utils.get_centroid: sum the raw counts over the pool,
%   normalise that one profile to 1e4 and log1p it, which weights each cell
%   by its depth rather than equally.
%
%   WHAT THIS COSTS IN ACCURACY. A settled pool reports one type for all of
%   its cells, so a type present as a few cells inside a pool is absorbed
%   rather than reported. A pseudobulk hides such a minority more thoroughly
%   than a sample would, because those cells' counts are averaged away
%   instead of showing up as a dissenting label. Rare-cell-type discovery is
%   the one job this mode is wrong for; use per-cell annotation.
%
%   INPUTS:
%     gid           - NumCells-by-1 pool ids, e.g. from
%                     PKG.I_PARTITIONCELLS
%     stage1        - one row per pool, with variables:
%                       Pool       - pool id as a string, matching GID
%                       Label      - the pool pseudobulk's label
%                       Confidence - that label's confidence in 0..1
%                       HalfLabel  - npools-by-2 string, the two half
%                                    pseudobulks' labels
%     annotatecells - function handle, [labels, weights] = annotatecells(idx).
%                     Called once, with a sorted column of cell indices,
%                     and only if some pool is unsettled. Must return one
%                     label per index, in order, and either one weight per
%                     index or [] for an unweighted vote.
%
%                     TWO OUTPUTS ARE REQUIRED, not optional. Probing for
%                     a second output by calling the handle and catching
%                     the failure does not work here: a named function
%                     with one output RUNS ITS BODY before MATLAB raises
%                     MATLAB:unassignedOutputs, so the fallback would
%                     re-run the classifier - minutes of Python, twice.
%                     A single-output handle therefore gets a clear error
%                     rather than a silent retry. Wrap it in
%                     @(idx) deal(yourfun(idx), []) if it has no weights.
%
%   NAME-VALUE:
%     MinConfidence      - (0.5) confidence a pool's pseudobulk label needs
%                          to settle the pool without labelling its cells
%     MinShare           - (0.6) share of a fully labelled pool its winning
%                          label must hold to settle it
%     RequireHalvesAgree - (true) also require both half pseudobulks to
%                          carry the whole pool's label
%
%   OUTPUTS:
%     labels  - NumCells-by-1 string. Settled pools carry their pool label;
%               mixed pools carry their cells' own labels.
%     report  - one row per pool: Pool, Size, Labelled, Label, Confidence,
%               HalvesAgree, Share, NDistinct, Status. Status is
%               "consensus" (settled on the pseudobulk), "escalated"
%               (cells labelled, then settled) or "mixed" (cells labelled,
%               not settled). Share and NDistinct are NaN for a pool whose
%               cells were never labelled - there is nothing to take a
%               share of, and reporting 1 there would read as certainty.
%     percell - NumCells-by-1 string of the labels the classifier returned
%               for individual cells, "" for the cells it was never asked
%               about. On a run where nothing escalated this is all "",
%               which is the point of the mode.
%     percellweight - NumCells-by-1 of the weights that came back with
%               those labels, NaN where none were asked for or none were
%               supplied. NaN rather than 1, so "no measurement" cannot be
%               mistaken for "full confidence".
%
% See also PKG.I_PARTITIONCELLS, PKG.I_MAJORITYVOTE, RUN.PY_SCIMILARITY.

arguments
    gid {mustBeNonempty}
    stage1 table
    annotatecells (1,1) function_handle
    opts.MinConfidence (1,1) double ...
        {mustBeGreaterThanOrEqual(opts.MinConfidence, 0), ...
        mustBeLessThanOrEqual(opts.MinConfidence, 1)} = 0.5
    opts.MinShare (1,1) double {mustBeGreaterThanOrEqual(opts.MinShare, 0), ...
        mustBeLessThanOrEqual(opts.MinShare, 1)} = 0.6
    opts.RequireHalvesAgree (1,1) logical = true
end

groups = string(gid(:));
ncells = numel(groups);
percell = strings(ncells, 1);
percellweight = nan(ncells, 1);

i_checkstage1(stage1, groups);

% Put STAGE1 in the pools' sorted order, so every vector below indexes the
% same way and the caller is free to hand its rows over in any order.
[poolname, ~, poolindex] = unique(groups);
[~, wherepool] = ismember(poolname, string(stage1.Pool));
stage1 = stage1(wherepool, :);
npools = numel(poolname);
poolsize = accumarray(poolindex, 1);

halvesagree = all(string(stage1.HalfLabel) == string(stage1.Label), 2);
settled = stage1.Confidence >= opts.MinConfidence;
if opts.RequireHalvesAgree
    settled = settled & halvesagree;
end

% Every cell of every unsettled pool, in one call.
if any(~settled)
    rest = find(ismember(groups, poolname(~settled)));
    [percell(rest), restweight] = i_ask(annotatecells, rest);
    if ~isempty(restweight)
        percellweight(rest) = restweight;
    end
end

report = table(poolname, poolsize, ...
    accumarray(poolindex, double(strlength(percell) > 0)), ...
    string(stage1.Label), stage1.Confidence, halvesagree, ...
    nan(npools, 1), nan(npools, 1), repmat("consensus", npools, 1), ...
    'VariableNames', {'Pool', 'Size', 'Labelled', 'Label', 'Confidence', ...
    'HalvesAgree', 'Share', 'NDistinct', 'Status'});

if any(~settled)
    labelled = strlength(percell) > 0;
    % One batched call produced every one of these labels, so the weights
    % are either all present or all absent; [] then means an unweighted
    % vote rather than a half-weighted one.
    voteweight = percellweight(labelled);
    if all(isnan(voteweight)), voteweight = []; end
    [~, ~, tally] = pkg.i_majorityvote(percell(labelled), ...
        groups(labelled), voteweight);
    [~, wheretally] = ismember(tally.Group, poolname);

    report.Label(wheretally) = tally.Label;
    report.Share(wheretally) = tally.Share;
    report.NDistinct(wheretally) = tally.NDistinct;
    report.Status(wheretally) = "escalated";
    report.Status(wheretally(tally.Share < opts.MinShare)) = "mixed";
end

labels = report.Label(poolindex);

% A pool the full labelling did not settle is not flattened: its cells keep
% what the classifier said about each of them. They all have a label,
% because reaching "mixed" means the pool was labelled in full.
mixedcells = ismember(groups, report.Pool(report.Status == "mixed"));
labels(mixedcells) = percell(mixedcells);
end


function i_checkstage1(stage1, groups)
% STAGE1 must describe exactly the pools GID contains. A mismatch would
% quietly label cells from whichever row ISMEMBER happened to land on, so
% it is an error rather than a warning.

needed = {'Pool', 'Label', 'Confidence', 'HalfLabel'};
missing = needed(~ismember(needed, stage1.Properties.VariableNames));
if ~isempty(missing)
    error('pkg:i_pooledannotate:stage1Variables', ...
        'STAGE1 is missing the variable(s) %s.', strjoin(missing, ', '));
end

if size(stage1.HalfLabel, 2) ~= 2
    error('pkg:i_pooledannotate:halfLabelWidth', ...
        ['STAGE1.HalfLabel must have 2 columns, one per half pool, but ', ...
        'has %d.'], size(stage1.HalfLabel, 2));
end

pools = unique(groups);
listed = string(stage1.Pool);
if numel(listed) ~= numel(unique(listed))
    error('pkg:i_pooledannotate:duplicatePool', ...
        'STAGE1.Pool lists the same pool more than once.');
end
if ~isequal(sort(listed(:)), sort(pools(:)))
    error('pkg:i_pooledannotate:poolMismatch', ...
        ['STAGE1 describes %d pool(s) and GID contains %d. They must be ', ...
        'the same set of pools.'], numel(listed), numel(pools));
end
end


function [labels, weights] = i_ask(annotatecells, idx)
% One batched call, with the contract checked. A classifier that silently
% returned a different number of labels than it was asked for would shift
% every label onto the wrong cell.
%
% The failure below is reported, never retried. See the ANNOTATECELLS note
% in the help: a one-output named function has already run its body by the
% time MATLAB complains, so calling it again would repeat the work.

try
    [labels, weights] = annotatecells(idx);
catch ME
    if ismember(ME.identifier, {'MATLAB:maxlhs', ...
            'MATLAB:unassignedOutputs', 'MATLAB:TooManyOutputs'})
        error('pkg:i_pooledannotate:annotatorOutputs', ...
            ['ANNOTATECELLS must return [labels, weights]. It returned ', ...
            'one output. Wrap it as @(idx) deal(yourfun(idx), []) if it ', ...
            'has no confidence to report.']);
    end
    rethrow(ME);
end

labels = reshape(string(labels), [], 1);
if numel(labels) ~= numel(idx)
    error('pkg:i_pooledannotate:labelCount', ...
        ['The annotate function was asked about %d cells and returned ', ...
        '%d labels. It must return one label per index, in order.'], ...
        numel(idx), numel(labels));
end
labels(ismissing(labels)) = "";

if ~isempty(weights)
    weights = reshape(double(weights), [], 1);
    if numel(weights) ~= numel(idx)
        error('pkg:i_pooledannotate:weightCount', ...
            ['The annotate function was asked about %d cells and ', ...
            'returned %d weights. Return one per index, or [].'], ...
            numel(idx), numel(weights));
    end
end
end
