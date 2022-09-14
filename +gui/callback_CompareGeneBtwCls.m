function callback_CompareGeneBtwCls(src,~)
    answer = questdlg('This function generate violinplot to show differences between cell groups. Continue?','');
    if ~strcmp(answer,'Yes'), return; end    

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc]=gui.i_select1class(sce);
    if isempty(thisc), return; end


    
    
selitems={'Expression of Gene','Library Size', ...
    'Predefined Cell Score','Other Attribute'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',selitems);
if tf1~=1, return; end

%try
    switch selitems{indx1}
        case 'Library Size'
            y=sum(sce.X);
            ttxt='Library Size';
        case 'Expression of Gene'
            [glist]=gui.i_selectngenes(sce);
            if isempty(glist)
                helpdlg('No gene selected.','');
                return;
            end
            [Xt]=gui.i_transformx(sce.X);
            % [cL]=i_getgrouporder(thisc);
            [~,cL,noanswer]=gui.i_reordergroups(thisc);
            if noanswer, return; end
            
            [answerc]=questdlg('Color violin plot?','');
            switch answerc
                case 'Yes'
                    colorit=true;
                case 'No'
                    colorit=false;
                otherwise
                    return;
            end

            gui.i_cascadeviolin(sce,Xt,thisc,glist, ...
                selitems{indx1},cL,colorit);


%         figure;
%         [~,idx]=ismember(glist,sce.g);
%         h=heatmap(Xt(idx,1:500));
%         h.Title = 'Gene Heatmap';
%         h.XLabel = 'Group';
%         h.YLabel = 'Genes';
%         h.Colormap = parula;
%         h.GridVisible = 'off';

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
    
    [~,cL,noanswer]=gui.i_reordergroups(thisc);
    if noanswer, return; end
    %[cL]=i_getgrouporder(thisc);
    f = figure('visible','off');
    pkg.i_violinplot(y,thisc,false,cL);
    title(strrep(ttxt,'_','\_'));
    ylabel(selitems{indx1});
    movegui(f,'center');
    set(f,'visible','on');
%catch ME
%    errordlg(ME.message);
%end

end



% function [cL]=i_getgrouporder(thisc)
%     [c,cL]=grp2idx(thisc);
%     [answer]=questdlg('Manually order groups?','');
%     switch answer
%         case 'Yes'
%             [newidx]=gui.i_selmultidlg(cL,sort(cL));
%             if length(newidx)~=length(cL)
%                 waitfor(helpdlg('Reordering is canceled or incomplete. Default order is used. Click OK to continue.',''));
%                 return;
%             end
%             cx=c;
%             for k=1:length(newidx)
%                 c(cx==newidx(k))=k;
%             end
%             cL=cL(newidx);
%         otherwise
%     end
% end
