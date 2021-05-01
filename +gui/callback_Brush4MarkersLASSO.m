function callback_Brush4MarkersLASSO(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    assert(isequal(FigureHandle.Children, FigureHandle.findobj('type','Axes')))
    
    axesh=FigureHandle.Children(1);    
    [ax,bx]=view(axesh);    
    assert(isequal(axesh.findobj('type','Scatter'),...
        FigureHandle.findobj('type','Scatter')))    
    h=axesh.Children(1);
    ptsSelected = logical(h.BrushData.');    
    
    if ~any(ptsSelected)
        warndlg("No cells are selected.");
        return;
    end
    assignin('base','ptsSelected',ptsSelected);
        
    [c,cL]=grp2idx(sce.c);
    if isscalar(unique(c))
        methodtag=1;
    else
        answer = questdlg('Select brushed cell group?');
        switch answer
            case 'Yes'
                uptsSelected=unique(c(ptsSelected));
                if isscalar(uptsSelected)
                    methodtag=2;   % whole group
                    ptsSelected=c==uptsSelected;
                else
                    errordlg('More than one group of brushed cells');
                    return;
                end
            case 'No'
                methodtag=1;       % only selected cells 
            otherwise
                return;
        end
    end    
    [numfig]=gui.gui_inputdlg;
    fw=gui.gui_waitbar;
    y=double(ptsSelected);
    sce.c=1+ptsSelected;
    X=sce.X';
    if issparse(X) 
        X=full(X); 
    end
    [B]=lasso(X,y,'DFmax',numfig*3,'MaxIter',1e3);
    % assignin('base','A',A);
    [~,ix]=min(abs(sum(B>0)-numfig));
    b=B(:,ix);
    idx=b>0;
    gui.gui_waitbar(fw);
    
    if ~any(idx)
        errordlg('No marker gene found')
        return;
    else
        fw=gui.gui_waitbar;
        markerlist=sce.g(idx);
        [~,jx]=sort(b(idx),'descend');
        markerlist=markerlist(jx);
        

        
        htmlfilename=cL{unique(c(ptsSelected))};
        pkg.i_markergeneshtml(sce,markerlist,numfig,...
                   [ax bx],htmlfilename,ptsSelected);
        gui.gui_waitbar(fw);        
    end
        fprintf('%d marker genes: ',length(markerlist));
        fprintf('%s ',markerlist)
        fprintf('\n')
%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end
