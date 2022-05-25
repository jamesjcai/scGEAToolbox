function callback_CompareGeneBtwCls(src,~)
    answer = questdlg('This function generate violinplot to show differences between cell groups. Continue?','');
    if ~strcmp(answer,'Yes'), return; end    

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc]=gui.i_select1class(sce);
    if isempty(thisc), return; end
    
a={'Expression of Selected Genes','Library Size','Predefined Cell Score','Other Attribute'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',a);
if tf1~=1, return; end

try
    switch a{indx1}
        case 'Library Size'
            y=sum(sce.X);
            ttxt='Library Size';
        case 'Expression of Selected Genes'
            [glist]=gui.i_selectngenes(sce);
            if isempty(glist)
                helpdlg('No gene selected.','');
                return;
            end
            [Xt]=gui.i_transformx(sce.X);        
            for k=1:length(glist)
                [~,idx]=ismember(glist(k),sce.g);
                y=full(Xt(idx,:));
                ttxt=sce.g(idx);
                
                f = figure('visible','off');
                pkg.i_violinplot(y,thisc);
                title(strrep(ttxt,'_','\_'));
                ylabel(a{indx1});
                tb=uitoolbar(f);
                pkg.i_addbutton2fig(tb,'off',{@i_savedata,y,thisc},'export.gif','Export data...');
                P = get(f,'Position');
                set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
                set(f,'visible','on');                
                drawnow;
            end
            return;
        case 'Predefined Cell Score'
            [~,T]=pkg.e_cellscores(sce.X,sce.g,0);
            listitems=T.ScoreType;
            [indx2,tf2] = listdlg('PromptString',...
                {'Select Class','',''},...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end
            fw=gui.gui_waitbar;
            [y]=pkg.e_cellscores(sce.X,sce.g,indx2);
            ttxt=T.ScoreType(indx2);
            gui.gui_waitbar(fw);
        case 'Other Attribute'
            [y,clable,~,newpickclable]=gui.i_select1state(sce,true);
            if isempty(y)
                helpdlg('No cell attribute is available.');
                return; 
            end
            if ~isempty(newpickclable)
                ttxt=newpickclable;
            else
                ttxt=clable;
            end
        otherwise
            return;
    end
    f = figure('visible','off');
    pkg.i_violinplot(y,thisc);
    title(strrep(ttxt,'_','\_'));
    ylabel(a{indx1});
    movegui(f,'center');
    set(f,'visible','on');
catch ME
    errordlg(ME.message);
end

end


function i_savedata(~,~,a,b)
    T=table(a(:),b(:));    
    T.Properties.VariableNames={'ExprLevel','GroupID'};
    T=sortrows(T,'ExprLevel','descend');
    T=sortrows(T,'GroupID');
    gui.i_exporttable(T,true);
end

