function [Tup, Tdn, paramset, usedset] = e_processdetable(T, paramset, parentfig)
%E_PROCESSDETABLE Filter a DE table into up- and down-regulated genes.
%
%   [TUP, TDN] = PKG.E_PROCESSDETABLE(T, PARAMSET, PARENTFIG)
%   [TUP, TDN, PARAMSET, USEDSET] = PKG.E_PROCESSDETABLE(...)
%
%   PARAMSET is {minDiffPct, minAbsLog2FC, adjPCutoff, sortBy, cutoffMode}
%   as returned by GUI.I_DEGPARAMSET. CUTOFFMODE is optional:
%       'Fixed' (default) - apply the three cutoffs as given.
%       'Automatic'       - if neither list reaches MINGENES (20) genes,
%                           relax the cutoffs one step at a time, never
%                           tightening any of them, until one list does or
%                           there is nothing left to relax. See
%                           IN_RELAXSTEPS for the order. If that is still
%                           not enough at adjusted P 0.05, the P cutoff is
%                           set exactly at the 20th gene of the better
%                           direction (IN_EXACTCUT), which adds a raw
%                           P-value cutoff to break BH ties.
%   An empty PARAMSET opens the dialog.
%
%   USEDSET holds the cutoffs actually applied, so callers can report them:
%   {minDiffPct, minAbsLog2FC, adjPCutoff, sortBy, cutoffMode, rawPCutoff},
%   where rawPCutoff is Inf unless the last automatic stage set it.

if nargin<3, parentfig = []; end
Tup=[];
Tdn=[];
usedset = [];
if nargin < 2, paramset = []; end
if isempty(paramset)
    [paramset] = gui.i_degparamset(false, parentfig);
end
if isempty(paramset)
   return;
else
    mindiffpct = paramset{1};
    minabsolfc = paramset{2};
    apvaluecut = paramset{3};
    sortbywhat = paramset{4};
end
cutoffmode = 'Fixed';
if numel(paramset) >= 5 && ~isempty(paramset{5})
    cutoffmode = paramset{5};
end

minGenes = 20;

a = T.(T.Properties.VariableNames{8}) - T.(T.Properties.VariableNames{7});

% Raw P-value cutoff; Inf means none. Only the last automatic stage sets it.
rawpcut = Inf;

[nup, ndn] = in_count(T, in_passing(T, a, mindiffpct, minabsolfc, apvaluecut, rawpcut));
if strcmp(cutoffmode, 'Automatic') && max(nup, ndn) < minGenes
    fprintf(['\nAutomatic cutoffs: %d up- and %d down-regulated genes pass ' ...
        'the entered cutoffs; relaxing until one list has >= %d.\n'], ...
        nup, ndn, minGenes);
    steps = in_relaxsteps();
    for k = 1:size(steps, 1)
        mindiffpct = min(mindiffpct, steps{k, 1});
        minabsolfc = min(minabsolfc, steps{k, 2});
        apvaluecut = max(apvaluecut, steps{k, 3});
        [nup, ndn] = in_count(T, in_passing(T, a, mindiffpct, minabsolfc, apvaluecut, rawpcut));
        fprintf('  step %d: diff(pct) >= %.2f, abs(log2FC) >= %.2f, adj. P <= %g -> %d up, %d down\n', ...
            k, mindiffpct, minabsolfc, apvaluecut, nup, ndn);
        if max(nup, ndn) >= minGenes
            break;
        end
    end

    if max(nup, ndn) < minGenes
        % Past significance. A ladder of adjusted-P values cannot stop near
        % MINGENES: BH ties most genes at one value when little is DE (on a
        % random split of the bundled data, 3,840 of 3,850 sat at 0.995),
        % so the next rung admits everything. BH is non-decreasing in the
        % raw P-value, so ranking on raw P is the adjusted-P ranking with
        % the ties broken. Cut at the MINGENES-th gene of whichever
        % direction gets there first, and apply that cutoff to both.
        [exactcut, rawpcut] = in_exactcut(T, a, mindiffpct, minabsolfc, minGenes);
        apvaluecut = max(apvaluecut, exactcut);
        [nup, ndn] = in_count(T, in_passing(T, a, mindiffpct, minabsolfc, apvaluecut, rawpcut));
        fprintf('  last step: adj. P <= %g and raw P <= %g -> %d up, %d down\n', ...
            apvaluecut, rawpcut, nup, ndn);
    end

    if max(nup, ndn) < minGenes
        warning('pkg:e_processdetable:tooFewGenes', ...
            ['Fewer than %d genes change in either direction even with ' ...
            'every cutoff relaxed; reporting the %d up and %d down that do.'], ...
            minGenes, nup, ndn);
    end
    if apvaluecut > 0.05
        warning('pkg:e_processdetable:notSignificant', ...
            ['Automatic cutoffs relaxed the adjusted P-value to %g. Genes ' ...
            'above 0.05 are not statistically significant; treat the lists ' ...
            'as a ranking, not as DE calls.'], apvaluecut);
    end
end
usedset = {mindiffpct, minabsolfc, apvaluecut, sortbywhat, cutoffmode, rawpcut};

isok = in_passing(T, a, mindiffpct, minabsolfc, apvaluecut, rawpcut);

fprintf(['\nDE genes with > %.2f%% difference in expression percentages, ' ...
'abs(log2FC) >= %.2f, and adjusted P-value < %.3f are retained.\n'], ...
mindiffpct*100,...
minabsolfc, apvaluecut);
if isfinite(rawpcut)
    fprintf('Raw P-value <= %g is also required (automatic cutoffs).\n', rawpcut);
end

Tup = T(T.avg_log2FC > 0 & isok, :);
Tdn = T(T.avg_log2FC < 0 & isok, :);

if ~isempty(sortbywhat)
    answer = sortbywhat;
else
    answer = gui.myQuestdlg(parentfig, 'Sort DE genes by adjusted P-value or fold change?','',...
        {'Adjusted P-value','Fold Change'},'Adjusted P-value');
end

switch answer
    case 'Adjusted P-value'
        Tup = sortrows(Tup, 'abs_log2FC', 'descend');
        Tup = sortrows(Tup, 'p_val_adj', 'ascend');
        Tdn = sortrows(Tdn, 'abs_log2FC', 'descend');
        Tdn = sortrows(Tdn, 'p_val_adj', 'ascend');
        disp('DE genes are sorted by adjusted P-value.');
    case 'Fold Change'
        Tup = sortrows(Tup, 'p_val_adj', 'ascend');
        Tup = sortrows(Tup, 'abs_log2FC', 'descend');
        Tdn = sortrows(Tdn, 'p_val_adj', 'ascend');
        Tdn = sortrows(Tdn, 'abs_log2FC', 'descend');
        disp('DE genes are sorted by absolute fold change (FC).');
    otherwise
        gui.myHelpdlg(parentfig, 'Keep DE gene tables unsorted.');
end
end

function isok = in_passing(T, a, mindiffpct, minabsolfc, apvaluecut, rawpcut)
isok = abs(a) >= mindiffpct & abs(T.avg_log2FC) >= minabsolfc & ...
       T.p_val_adj <= apvaluecut & T.p_val <= rawpcut;
end

function [apvaluecut, rawpcut] = in_exactcut(T, a, mindiffpct, minabsolfc, minGenes)
% The smallest raw-P cutoff at which one direction has MINGENES genes
% passing the effect-size cutoffs, and the largest adjusted P that admits.
% With fewer than MINGENES candidates in both directions, admit them all.
cand = in_passing(T, a, mindiffpct, minabsolfc, Inf, Inf);
pup = sort(T.p_val(cand & T.avg_log2FC > 0));
pdn = sort(T.p_val(cand & T.avg_log2FC < 0));
kth = Inf(1, 2);
if numel(pup) >= minGenes, kth(1) = pup(minGenes); end
if numel(pdn) >= minGenes, kth(2) = pdn(minGenes); end
rawpcut = min(kth);
apvaluecut = max([T.p_val_adj(cand & T.avg_log2FC ~= 0 & T.p_val <= rawpcut); 0]);
end

function [nup, ndn] = in_count(T, isok)
nup = nnz(isok & T.avg_log2FC > 0);
ndn = nnz(isok & T.avg_log2FC < 0);
end

function steps = in_relaxsteps()
% One row per step: {minDiffPct, minAbsLog2FC, adjPCutoff}. Each step only
% loosens: E_PROCESSDETABLE takes MIN of the first two and MAX of the third
% against what is already in force, so a cutoff the user already set looser
% is never tightened. Effect size goes first and the P-value last; nothing
% here goes past 0.05. Beyond that, IN_EXACTCUT picks the cutoff directly.
steps = {
    0.05, 0.585, 0.01    % 1.5-fold
    0.05, 0.585, 0.05
    0.05, 0.26,  0.05    % 1.2-fold
    0.01, 0.26,  0.05
    0,    0,     0.05    % significance alone
    };
end
