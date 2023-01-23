function callback_CompareGeneBtwCls(src,~)
    answer = questdlg('This function generate violinplot to show differences between cell groups. Continue?','');
    if ~strcmp(answer,'Yes'), return; end    

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc]=gui.i_select1class(sce);
    if isempty(thisc), return; end

    if length(unique(thisc))==1
        answer=questdlg("All cells are in the same group. Continue?","");
        switch answer
            case 'Yes'
            otherwise
                return;
        end
    end

selitems={'Expression of Gene', ...
    'Predefined Cell Score','TF Activity Score','Differentiation Potency',...
    '--------------------------------',...
    'Library Size','Other Attribute'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',selitems);
if tf1~=1, return; end

%try
    switch selitems{indx1}
        case 'Differentiation Potency'
            [a]=contains(sce.list_cell_attributes(1:2:end),'cell_potency');
            if ~any(a)
                answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                [yes,specisid]=ismember(lower(answer2),{'human','mouse'});
                if ~yes, return; end                
                sce=sce.estimatepotency(specisid);
            end
            [yes,idx]=ismember({'cell_potency'},sce.list_cell_attributes(1:2:end));
            if yes
                y=sce.list_cell_attributes{idx+1};
                ttxt='Differentiation Potency';
            else                
                return;
            end           
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
            
            colorit=true;
%             [answerc]=questdlg('Color violin plot?','');
%             switch answerc
%                 case 'Yes'
%                     colorit=true;
%                 case 'No'
%                     colorit=false;
%                 otherwise
%                     return;
%             end
            cL=strrep(cL,'_','\_');
            thisc=strrep(thisc,'_','\_');
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
            [indx2,tf2] = listdlg('PromptString','Select Class',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end
            fw=gui.gui_waitbar;
            y=pkg.e_cellscores(sce.X,sce.g,indx2);
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
        case 'TF Activity Score'
            [~,T]=pkg.e_tfactivityscores(sce.X,sce.g,0);
            listitems=unique(T.tf);
            [indx2,tf2] = listdlg('PromptString','Select Class',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end
            fw=gui.gui_waitbar;
            [cs,tflist]=sc_tfactivity(sce.X,sce.g,[]);
            idx=find(tflist==string(listitems{indx2}));
            assert(length(idx)==1)
            [y]=cs(idx,:);
            ttxt=listitems{indx2};
            gui.gui_waitbar(fw);
        otherwise
            return;
    end
    
    [~,cL,noanswer]=gui.i_reordergroups(thisc);
    if noanswer, return; end

%     [answerc]=questdlg('Color violin plot?','');
%     switch answerc
%         case 'Yes'
%             colorit=true;
%         case 'No'
%             colorit=false;
%         otherwise
%             return;
%     end

    colorit=true;
    %[cL]=i_getgrouporder(thisc);
    f=figure('visible','off');
    tb=uitoolbar(f);
    pkg.i_addbutton2fig(tb,'off',{@i_savedata,y,thisc}, ...
        'export.gif','Export data...');
    pkg.i_addbutton2fig(tb,'off',{@i_testdata,y,thisc}, ...
        'exportx.gif','ANOVA/T-test...');    
    pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3}, ...
        "powerpoint.gif",'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb,'off',@i_invertcolor, ...
        "xpowerpoint.gif",'Switch BW/Color');

    cL=strrep(cL,'_','\_');
    thisc=strrep(thisc,'_','\_');
    pkg.i_violinplot(y,thisc,colorit,cL);
    title(strrep(ttxt,'_','\_'));
    ylabel(selitems{indx1});
    movegui(f,'center');
    set(f,'visible','on');
    
%catch ME
%    errordlg(ME.message);
%end

    function i_invertcolor(~,~)
        colorit=~colorit;
        delete(gca);
        pkg.i_violinplot(y,thisc,colorit,cL);
    end

end

function i_savedata(~,~,a,b)
    T=table(a(:),b(:));    
    T.Properties.VariableNames={'ExprLevel','GroupID'};
    T=sortrows(T,'ExprLevel','descend');
    T=sortrows(T,'GroupID');
    gui.i_exporttable(T,true);
end

function i_testdata(~,~,y,grp)
    if length(unique(grp))==2
        id=grp2idx(grp);
        [~,p,~,stats] = ttest2(y(id==1),y(id==2));
        stats.p=p;
        tbl=struct2table(stats);        
    else
        [~,tbl] = anova1(y,grp,"on");
        tbl=cell2table(tbl);
    end
    gui.i_exporttable(tbl,true);
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
