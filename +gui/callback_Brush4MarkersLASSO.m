function callback_Brush4MarkersLASSO(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
%   assert(isequal(FigureHandle.Children,...
%         FigureHandle.findobj('type','Axes')))
    
    % axesh=FigureHandle.Children(1)
    axesh=FigureHandle.findobj('type','Axes');
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
    [ptsSelected,letdoit]=gui.i_expandbrushed(ptsSelected,sce);
    if ~letdoit, return; end
%     [c,cL]=grp2idx(sce.c);
%     if isscalar(unique(c))
%         % methodtag=1;
%     else
%         answer = questdlg(sprintf('Select brushed cells'' group?\nYES to select brushed cells'' group\nNO to select brushed cells only'));
%         switch answer
%             case 'Yes'
%                 uptsSelected=unique(c(ptsSelected));
%                 if isscalar(uptsSelected)
%                     % methodtag=2;   % whole group
%                     ptsSelected=c==uptsSelected;
%                 else
%                     errordlg('More than one group of brushed cells');
%                     return;
%                 end
%             case 'No'
%                 % methodtag=1;       % only selected cells 
%             otherwise
%                 return;
%         end
%     end


    [numfig]=gui.i_inputnumg;
    if isempty(numfig), return; end
    fw=gui.gui_waitbar;
    y=double(ptsSelected);
    sce.c=1+ptsSelected;
    X=sce.X';
    try
        if issparse(X) 
            X=full(X); 
        end
        [B]=lasso(X,y,'DFmax',numfig*3,'MaxIter',1e3);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        rethrow(ME);
    end
    
    [~,ix]=min(abs(sum(B>0)-numfig));
    b=B(:,ix);
    idx=b>0;
    gui.gui_waitbar(fw);
    
    if ~any(idx)
        errordlg('No marker gene found')
        return;
    else
         markerlist=sce.g(idx);
         [~,jx]=sort(b(idx),'descend');
         markerlist=markerlist(jx);
%        htmlfilename=cL{unique(c(ptsSelected))};
%        pkg.i_markergeneshtml(sce,markerlist,numfig,...
%                    [ax bx],htmlfilename,ptsSelected);        
        for kk=1:length(markerlist)
            gui.i_cascadefig(sce,markerlist(end-(kk-1)),ax,bx,kk);
        end       
    end
        fprintf('%d marker genes: ',length(markerlist));
        fprintf('%s ',markerlist)
        fprintf('\n')
%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end
