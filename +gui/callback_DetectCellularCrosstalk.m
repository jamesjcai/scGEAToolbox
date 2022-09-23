function callback_DetectCellularCrosstalk(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc,~]=gui.i_select1class(sce,false);
    if isempty(thisc), return; end
    
    [c,cL]=grp2idx(thisc);

    [idx]=gui.i_selmultidlg(cL);
    if isempty(idx), return; end
    if numel(idx)<2
        warndlg('Need at least 2 cell types');
        return;
    end
    selected=ismember(c,idx);
    fw=gui.gui_waitbar;
    sce=sce.selectcells(selected);
    

    [OUT]=run.talklr(sce.X,sce.g,cL(c(selected)));
    gui.gui_waitbar(fw);    
    
    n=length(OUT.ligandok);
    if n==0
        warndlg('Not detected.');
    end
    
    labels = {'Save OUT to variable named:'};
    vars = {'OUT'};
    values = {OUT};

    if ~(ismcc || isdeployed)
        [f,ft]=export2wsdlg(labels, vars, values);
        waitfor(f);
        if ft
            disp('Run >> gui.i_crosstalkgraph(OUT,k,sce); to plot crosstalk graph for ligand-receptor pair k.')
        end
    end

    listitems=cell(n,1);
    for k=1:n
        listitems{k}=sprintf('%s -> %s (KL = %.2f)',...
            OUT.ligandok(k), OUT.receptorok(k),...
            OUT.KL(k));
    end
    i_displyres(listitems);
    
function i_displyres(listitems)
    [indx2,tf2] = listdlg('PromptString',...
        {'Select ligand-receptor pairs to plot'},...
         'SelectionMode','single','ListString',listitems,...
         'ListSize',[210,300]);
     if tf2==1         
            kk=indx2;
            gui.i_crosstalkgraph(OUT,kk,sce);
            i_displyres(listitems);
     elseif tf2==0

     end
end



end