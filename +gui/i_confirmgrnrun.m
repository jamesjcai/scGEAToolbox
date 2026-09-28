function ok = i_confirmgrnrun(method, x, numNetworks, parentfig)
%I_CONFIRMGRNRUN Warn before a network run that is likely to disappoint.
%   OK = gui.i_confirmgrnrun(METHOD, X, NUMNETWORKS, parentfig) shows each
%   warning gui.i_grnrunwarnings finds for METHOD (a row of
%   net.grnmethods()) on the transformed genes-by-cells matrix X, and
%   returns false as soon as the user declines one.
%
%   See also gui.i_grnrunwarnings, gui.i_selectgrnmethod.

if nargin < 4, parentfig = []; end

ok = true;
for w = gui.i_grnrunwarnings(method, x, numNetworks)
    answer = gui.myQuestdlg(parentfig, char(w.Message), char(w.Title), ...
        {'Continue', 'Cancel'}, char(w.Default));
    ok = strcmp(answer, 'Continue');
    if ~ok, return; end
end
end
