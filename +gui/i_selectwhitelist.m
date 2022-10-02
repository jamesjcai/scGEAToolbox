function [whitelist]=i_selectwhitelist(sce)
           
    answer = questdlg('Genes in whitelist will not be removed. Select whitelist genes?',...
                'Whitelist Genes','Yes','No','Cancel','No');
    switch answer
        case 'Yes'
            whitelist=0;
%             [gsorted]=gui.i_sortgenenames(sce);
%             if isempty(gsorted), return; end
%             [idx]=gui.i_selmultidlg(gsorted);
%             if isempty(idx), return; end
%             if isscalar(idx) && idx==0, return; end
% whitelist=gsorted(idx);
            [whitelist]=gui.i_selectngenes(sce);
        case 'No'
            whitelist=[];     % when isempty, continue..
            return;
        case 'Cancel'
            whitelist=0;
            return;
        otherwise
            whitelist=0;
            return;
    end
end