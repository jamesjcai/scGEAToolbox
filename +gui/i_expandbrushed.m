function [ptsSelected, updated] = ...
    i_expandbrushed(ptsSelected, sce, FigureHandle)

if nargin < 3, FigureHandle = []; end
updated = false;

[c, ~] = grp2idx(sce.c);
if ~isscalar(unique(c)) && isscalar(unique(c(ptsSelected)))
    answer = gui.myQuestdlg(FigureHandle, sprintf('Select brushed cells only or expand to cell group?'), '', ...
        {'Brushed cells only', 'Expand to cell group'}, 'Brushed cells only');
    switch answer
        case 'Expand to cell group'
            uptsSelected = unique(c(ptsSelected));
            if isscalar(uptsSelected)
                % methodtag=2;   % whole group
                ptsSelected = c == uptsSelected;
                updated = true;
            else
                errordlg('More than one group of brushed cells');
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
