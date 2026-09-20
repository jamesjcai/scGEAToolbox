function [ptsSelected, updated] = ...
i_expandbrushed(ptsSelected, sce, parentfig)

if nargin < 3, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
updated = false;

[c, ~] = findgroups(sce.c);
if ~isscalar(unique(c)) && isscalar(unique(c(ptsSelected)))
    answer = gui.myQuestdlg(parentfig, sprintf('Select brushed cells only or expand to cell group?'), '', ...
        {'Brushed cells only', 'Expand to cell group'}, 'Brushed cells only');
    if isempty(answer), return; end
    switch answer
        case 'Expand to cell group'
            uptsSelected = unique(c(ptsSelected));
            if isscalar(uptsSelected)
                % methodtag=2;   % whole group
                ptsSelected = c == uptsSelected;
                updated = true;
                % Grow the highlight to match. The user was asked which
                % cells to work on and said the group; a plot still
                % showing the original handful says the answer was
                % ignored, and the handlers that go on to run for several
                % seconds give no other sign of what they are running on.
                % DELETE SELECTED CELLS only looked right here because the
                % cells it took then disappeared.
                in_showexpanded(parentfig, ptsSelected);
            else
                gui.myErrordlg(parentfig, 'More than one group of brushed cells');
                return;
            end
        case 'Brushed cells only'
            updated = true;
            % methodtag=1;       % only selected cells
        otherwise
            return;
    end
else
    updated = true;
end

end

function in_showexpanded(parentfig, ptsSelected)
% The plot's own scatter, when there is one to find and it is a plot of these
% cells. GUI.I_SETBRUSHDATA declines anything else.

if isempty(parentfig) || ~pkg.i_isvalid(parentfig), return; end
h = findobj(parentfig, 'Type', 'Scatter');
if isempty(h), return; end
if gui.i_setbrushdata(h(1), ptsSelected)
    % The caller's next move is usually a dialog or a long analysis, and
    % neither yields to the event queue on its own.
    drawnow;
end
end
