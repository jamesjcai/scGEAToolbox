function [numPermutations, fileTag] = i_dvpermutations(parentfig, method)
%I_DVPERMUTATIONS Ask whether SC_DVG's p-values come from permutations.
%
%   [NUMPERMUTATIONS, FILETAG] = GUI.I_DVPERMUTATIONS(PARENTFIG, METHOD)
%   asks the user and returns 0 (the closed-form p-value, SC_DVG's default)
%   or 30 permutations, or [] if the dialog was dismissed. FILETAG is '' or '_perm<K>', a file-name suffix so the two
%   kinds of run do not overwrite each other.
%
%   METHOD 'brennecke' returns 0 without asking: SC_DVG does not support
%   permutations for it.
%
%   See "Calibration" in SC_DVG for why the choice matters: the default
%   p-value is not calibrated on real cells, and every caller of this
%   dialog thresholds it at BH q <= 0.05.

arguments
    parentfig = []
    method = 'splinefit'
end

defaultCount = 30;
numPermutations = 0;
fileTag = '';
if strcmpi(method, 'brennecke')
    return;
end

% The explanation goes in the message, not on the buttons; see
% GUI.I_DVDIRECTION.
optFast = 'Standard (fast)';
optPerm = sprintf('Permutation (%d)', defaultCount);
message = {'How should DV p-values be computed?'; ''; ...
    [optFast ' -- one null for all genes; on real cells it passes ' ...
    'many genes that do not differ']; ...
    [optPerm ' -- calibrated per gene; about ' ...
    sprintf('%dx', defaultCount + 1) ' slower']};
answer = gui.myQuestdlg(parentfig, message, 'DV P-values', ...
    {optFast, optPerm}, optFast);
switch answer
    case optFast
        numPermutations = 0;
    case optPerm
        numPermutations = defaultCount;
        fileTag = sprintf('_perm%d', defaultCount);
    otherwise
        % Dismissed.
        numPermutations = [];
end
end
