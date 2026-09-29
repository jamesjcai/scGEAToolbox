function [best, agreement, n] = i_votelabels(answers)
%I_VOTELABELS Majority vote over repeated LLM answers, with an agreement score.
%
%   [best, agreement, n] = LLM.I_VOTELABELS(answers)
%
%   ANSWERS is a string array of independent replies to the same question.
%   BEST is the most common one, AGREEMENT its share of the votes, and N the
%   number of usable answers.
%
%   Votes are compared case- and punctuation-insensitively, so "T cells",
%   "T-cells" and "t cells" count as one answer; the most common ORIGINAL
%   spelling is what gets returned, since that is what a reader wants to see.
%
%   AGREEMENT is the point of this function as much as BEST. A label chosen
%   3/3 (1.00) and one chosen 1/3 because every reply differed (0.33) are
%   both "the answer", and only the score distinguishes a stable call from a
%   coin toss. An empty input gives "Unknown" with agreement 0, which is the
%   honest reading of "the model never answered".
%
%   See also LLM.E_CELLTYPEANNO, LLM.I_CLEANLLMLABEL.

arguments
    answers string
end

answers = answers(:);
answers = answers(~ismissing(answers) & strlength(strtrim(answers)) > 0);
n = numel(answers);

if n == 0
    best = "Unknown";
    agreement = 0;
    return
end

% 'stable' orders the distinct answers by first appearance, so max() breaks
% a tie in favour of the one given FIRST rather than the one that happens to
% sort first alphabetically. Callers control that order - repeated samples
% arrive in sampling order, and sc_annotatecells passes methods in priority
% order - so a tie resolves to something explainable instead of arbitrary.
key = lower(regexprep(answers, '[^a-zA-Z0-9]', ''));
[~, ~, idx] = unique(key, 'stable');
counts = accumarray(idx, 1);
[topCount, winner] = max(counts);

% Among the replies that agree, return the spelling used most often.
members = answers(idx == winner);
[uMembers, ~, mIdx] = unique(members, 'stable');
mCounts = accumarray(mIdx, 1);
[~, mWinner] = max(mCounts);
best = uMembers(mWinner);

agreement = topCount / n;
end
