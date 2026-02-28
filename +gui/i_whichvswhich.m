function [i1, i2, cL1, cL2, cancelled] = i_whichvswhich(FigureHandle, i1, i2, cL1, cL2)
%I_WHICHVSWHICH Present "Which vs. which?" dialog and swap group order if requested.
%
%   [i1, i2, cL1, cL2, cancelled] = gui.i_whichvswhich(FigureHandle, i1, i2, cL1, cL2)
%
%   Returns cancelled=true if the user dismissed the dialog.

cancelled = false;
a = sprintf('%s vs. %s', cL1{1}, cL2{1});
b = sprintf('%s vs. %s', cL2{1}, cL1{1});
answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?', '', {a, b}, a);
switch answer
    case a
        % keep original order
    case b
        [i1, i2] = deal(i2, i1);
        [cL1, cL2] = deal(cL2, cL1);
    otherwise
        cancelled = true;
end
end
