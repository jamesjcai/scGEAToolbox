function [ptsSelected,updated]=i_expandbrushed(ptsSelected,sce)
updated=false;
        
    [c,~]=grp2idx(sce.c);
    if ~isscalar(unique(c)) && isscalar(unique(c(ptsSelected)))
        answer = questdlg(sprintf('Select brushed cells'' group?\nYES to select brushed cells'' group\nNO to select brushed cells only'));
        switch answer
            case 'Yes'
                uptsSelected=unique(c(ptsSelected));
                if isscalar(uptsSelected)
                    % methodtag=2;   % whole group
                    ptsSelected=c==uptsSelected;
                    updated=true;
                else
                    errordlg('More than one group of brushed cells');
                    return;
                end
            case 'No'
                updated=true;
                % methodtag=1;       % only selected cells 
            otherwise
                return;
        end
    else
        updated=true;       
    end

end
