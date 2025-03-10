function [whitelist] = i_selectwhitelist(sce, parentfig)

if nargin<2, parentfig=[]; end

answer = gui.myQuestdlg(parentfig, 'Genes in whitelist will not be removed. Select whitelist genes?', ...
    'Whitelist Genes', {'Yes', 'No', 'Cancel'}, 'Yes');
switch answer
    case 'Yes'
        %            whitelist=0;
        %             [gsorted]=gui.i_sortgenenames(sce);
        %             if isempty(gsorted), return; end
        %             [idx]=gui.i_selmultidlg(gsorted);
        %             if isempty(idx), return; end
        %             if isscalar(idx) && idx==0, return; end
        % whitelist=gsorted(idx);
        [whitelist] = gui.i_selectngenes(sce, [], parentfig);
        if isempty(whitelist)
            waitfor(gui.myHelpdlg(parentfig, 'No whitelist gene selected. Click OK to continue.', ''));
            return;
        end
        [y] = ismember(whitelist, sce.g);
        if ~all(y)
            whitelist = whitelist(y);
        end
        if ~isempty(whitelist)
            waitfor(gui.myHelpdlg(parentfig, sprintf('%d whitelist gene(s) selected. Click OK to continue.', ...
                length(whitelist)), ''));
        else
            waitfor(gui.myHelpdlg(parentfig, 'No whitelist gene selected. Click OK to continue.', ''));
        end
    case 'No'
        whitelist = []; % when isempty, continue..
        return;
    case 'Cancel'
        whitelist = 0;
        return;
    otherwise
        whitelist = 0;
        return;
end
end