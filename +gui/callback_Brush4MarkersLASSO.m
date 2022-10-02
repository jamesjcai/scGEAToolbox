function callback_Brush4MarkersLASSO(src,~,sce)
        FigureHandle=src.Parent.Parent;

    if nargin<3
        sce=guidata(FigureHandle);
    end

%    FigureHandle=src.Parent.Parent;
%    sce=guidata(FigureHandle);
%   assert(isequal(FigureHandle.Children,...
%         FigureHandle.findobj('type','Axes')))
    
    % axesh=FigureHandle.Children(1)
    axesh=FigureHandle.findobj('type','Axes');
    [ax,bx]=view(axesh);    
    assert(isequal(axesh.findobj('type','Scatter'),...
        FigureHandle.findobj('type','Scatter')))    
    %axesh.Children(1)
    %isequal(axesh.findobj('type','Scatter'),axesh.Children(2))
    h=axesh.findobj('type','Scatter');
    ptsSelected = logical(h.BrushData.');
    

    if ~any(ptsSelected)
        %warndlg("No cells are selected.");
        %return;
        [ptsSelected]=gui.i_select1classcells(sce,false);
        if isempty(ptsSelected), return; end
    else
        assignin('base','ptsSelected',ptsSelected);
        [ptsSelected,letdoit]=gui.i_expandbrushed(ptsSelected,sce);
        if ~letdoit, return; end
    end


    [numfig]=gui.i_inputnumg(500);
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
        errordlg('No marker found')
        return;
    else
         markerlist=sce.g(idx);
         [~,jx]=sort(b(idx),'descend');
         markerlist=markerlist(jx);

        fprintf('%d marker genes: ',length(markerlist));
        fprintf('%s ',markerlist)
        fprintf('\n')

        gui.i_exporttable(table(markerlist),true, ...
            'T',"markerlist");


        [answer]=questdlg('Plot expression of markers?');
        if isempty(answer), return; end
        switch answer
            case 'Yes'
                [methodid]=gui.i_pickscatterstem('Scatter');
                if isempty(methodid), return; end
                F=cell(length(markerlist),1);
                for kk=1:length(markerlist)
                    F{kk}=gui.i_cascadefig(sce,markerlist(end-(kk-1)), ...
                        ax,bx,kk,methodid);            
                end
                gui.i_export2pptx(F,flipud(markerlist(:)));                
        end
    end
%     pause(2);
%     export2wsdlg({'Save marker list to variable named:'},...
%             {'g_markerlist'},{markerlist});
end
