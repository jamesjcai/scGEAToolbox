function txt = i_stashnotice(stashname)
%I_STASHNOTICE Sentence saying where the replaced cell type labels were kept.
% STASHNAME is what PKG.I_STASHCELLTYPEHISTORY returned. Returns "" when it
% stashed nothing, so a caller can append this to whatever it already
% reports without testing first:
%
%   gui.myHelpdlg(fig, msg + gui.i_stashnotice(stashname));
%
% The text opens with a blank line, so it reads as its own paragraph after
% the caller's own summary.
%
% Each route is given as a full menu path because the three live under two
% different menus - Compare and Assign under Annotate, View/Export under Edit
% - so naming the items alone sent the reader hunting in the wrong menu.
%
% The pre-flight warning (GUI.I_CONFIRMOVERWRITECELLTYPE) has to name the
% attribute generically, because the number is not known until the stash
% happens. This names the real one, after the fact, so the labels can
% actually be found again.
%
% See also PKG.I_STASHCELLTYPEHISTORY, GUI.I_CONFIRMOVERWRITECELLTYPE.

txt = "";
if nargin < 1 || isempty(stashname) || strlength(string(stashname)) == 0
    return;
end
txt = string(sprintf(['\n\nThe previous cell type labels were kept ', ...
    'as the ''%s'' cell attribute. Annotate > Compare Cell Type ', ...
    'Annotations shows it next to the current one. Annotate > Assign ', ...
    'Cell Type from Cell Attribute Table switches back to it. Edit > ', ...
    'View/Export Cell Attribute Table lists the values.'], ...
    string(stashname)));
end
