function [voted, nchanged, tally] = i_majorityvote(labels, groups, weights)
%I_MAJORITYVOTE Replace per-item labels with the most common label per group.
%
%   voted = pkg.i_majorityvote(labels, groups)
%   voted = pkg.i_majorityvote(labels, groups, weights)
%   [voted, nchanged, tally] = pkg.i_majorityvote(___)
%
%   Reference-model annotation (GUI.CALLBACK_RUNSCIMILARITY,
%   GUI.CALLBACK_RUNPANHUMANPY) embeds and classifies every cell
%   independently, so a cluster that is plainly one cell type still comes
%   back speckled with a few percent of its neighbours' types. This
%   collapses that: every cell in a group takes the label that most of the
%   group's cells were given.
%
%   Empty and missing labels take no part in the vote and are never
%   returned as a winner, so a handful of unlabelled cells cannot decide a
%   group and a group with nothing to count keeps what it had. Ties go to
%   the alphabetically first of the tied labels, so the same input always
%   gives the same answer.
%
%   WEIGHTS makes it a weighted vote: each item contributes its weight
%   rather than one vote, and the winner is the label holding the most
%   weight. The point is that a per-cell classifier says how sure it was -
%   SCimilarity returns vsAll_weighted, the winning label's share of the 50
%   reference neighbours - and an unweighted vote throws that away, letting
%   eight barely-decided cells outvote seven certain ones. An item with
%   zero weight cannot win and cannot help a label win; a group in which no
%   candidate carries positive weight is left alone, exactly like a group
%   with no usable label.
%
%   Weighting changes which label wins. It does not turn a mixed group into
%   a pure one, and a group that is genuinely half one type and half
%   another still gets flattened onto whichever half was more confident -
%   read SHARE before believing the answer.
%
%   The vote is only as good as the grouping. Passing an under-clustered
%   partition merges real cell types into whichever one is larger, and
%   that loss is not recoverable from VOTED - keep the per-item labels.
%
%   INPUTS:
%     labels  - N-by-1 labels (string, cellstr, categorical or numeric)
%     groups  - N-by-1 group ids of the same length, e.g. SCE.C_CLUSTER_ID
%     weights - N-by-1 non-negative weights, or [] for an unweighted vote.
%               [] and a vector of ones give identical results, down to
%               the tie-breaking.
%
%   OUTPUTS:
%     voted    - N-by-1 string, LABELS with each group collapsed to its
%                winner
%     nchanged - number of entries in which VOTED differs from LABELS
%     tally    - one row per group: the group id, its winning label, how
%                many of its members carried that label, the group size,
%                the winner's share, and how many distinct labels the
%                group carried. A share near 0.5 is a group the vote
%                barely decided.
%
%                COUNT counts items and SHARE measures weight, so under a
%                weighted vote SHARE is not COUNT/SIZE. A label carried by
%                5 confident cells can beat one carried by 10 doubtful
%                ones, and the row then reads Count 5, Size 15, Share
%                0.67. The two coincide exactly when WEIGHTS is omitted.
%
%                NDISTINCT and SHARE are not interchangeable as a measure
%                of how mixed a group is, and which one to read depends on
%                how many of its members were labelled. Three distinct
%                labels among 5 cells is surprising; three among 50 is
%                what per-cell classifier noise produces in a pure group.
%                Count distinct labels on a small sample, and take the
%                share on a large one.
%
% See also GUI.CALLBACK_RUNSCIMILARITY, PKG.I_STASHCELLTYPEHISTORY.

arguments
    labels
    groups
    weights {mustBeNumeric, mustBeReal} = []
end

labels = i_tocolumn(labels);
groups = i_tocolumn(groups);

if numel(labels) ~= numel(groups)
    error('pkg:i_majorityvote:sizeMismatch', ...
        'LABELS has %d entries and GROUPS has %d. Pass one group id per label.', ...
        numel(labels), numel(groups));
end

if isempty(weights)
    weights = ones(numel(labels), 1);
else
    weights = double(reshape(weights, [], 1));
    if numel(weights) ~= numel(labels)
        error('pkg:i_majorityvote:weightMismatch', ...
            'WEIGHTS has %d entries and LABELS has %d. Pass one weight per label.', ...
            numel(weights), numel(labels));
    end
    if any(~isfinite(weights)) || any(weights < 0)
        error('pkg:i_majorityvote:badWeight', ...
            ['WEIGHTS must be finite and non-negative. A NaN weight would ', ...
            'poison its whole group''s total, and a negative one would ', ...
            'let an item vote against a label.']);
    end
end

voted = labels;
if isempty(labels)
    nchanged = 0;
    tally = table(strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), 'VariableNames', ...
        {'Group', 'Label', 'Count', 'Size', 'Share', 'NDistinct'});
    return;
end

[gid, gname] = findgroups(groups);
ngroups = numel(gname);
winner = strings(ngroups, 1);
wincount = zeros(ngroups, 1);
groupsize = zeros(ngroups, 1);
ndistinct = zeros(ngroups, 1);
share = zeros(ngroups, 1);

for k = 1:ngroups
    rows = gid == k;
    groupsize(k) = sum(rows);
    % The denominator is the group's whole weight, unlabelled items
    % included, so a group half of which came back blank cannot report a
    % share above one half. That is the unweighted behaviour too, where
    % SHARE was the winning count over the group size rather than over the
    % number labelled.
    grouptotal = sum(weights(rows));

    keep = rows & strlength(labels) > 0 & weights > 0;
    candidates = labels(keep);
    if isempty(candidates) || grouptotal == 0, continue; end

    % UNIQUE sorts, so MAX picking the first maximum resolves a tie to the
    % alphabetically first label rather than to input order.
    [uniquelabels, ~, idx] = unique(candidates);
    counts = accumarray(idx, 1);
    weighted = accumarray(idx, weights(keep));
    ndistinct(k) = numel(uniquelabels);
    [winweight, w] = max(weighted);
    wincount(k) = counts(w);
    share(k) = winweight/grouptotal;
    winner(k) = uniquelabels(w);
    voted(rows) = winner(k);
end

nchanged = sum(voted ~= labels);
tally = table(gname(:), winner, wincount, groupsize, share, ...
    ndistinct, 'VariableNames', ...
    {'Group', 'Label', 'Count', 'Size', 'Share', 'NDistinct'});
end


function s = i_tocolumn(c)
% Any of the c_* attribute types as a column of strings, with missing
% values flattened to "" so they read as unlabelled rather than as a label
% of their own. Mirrors GUI.I_GETCELLGROUPS.

if isempty(c)
    s = strings(0, 1);
    return;
end
s = string(c);
s = reshape(s, [], 1);
s(ismissing(s)) = "";
s = strip(s);
end
