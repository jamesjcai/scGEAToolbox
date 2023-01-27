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
    'TF Activity Score','Differentiation Potency',...
    'MSigDB Signature Score',...    
    '--------------------------------',...
    'Predefined Cell Score',...
    'Define New Score...',...
    '--------------------------------',...
    'Library Size','Other Attribute'};
[indx1,tf1]=listdlg('PromptString',...
    'Select a metric for comparison.',...
    'SelectionMode','single','ListString',selitems);
if tf1~=1, return; end

%try
    switch selitems{indx1}
        case 'Define New Score...'
            ttxt='Customized Score';
            [posg]=gui.i_selectngenes(sce.g);
            if isempty(posg)
                helpdlg('No feature genes selected.','')
                return;
            end
            [y]=gui.e_cellscore(sce,posg);            
        case 'MSigDB Signature Score'
            stag=gui.i_selectspecies(2,true);
            if isempty(stag), return; end            
            try
                [posg,ctselected]=gui.i_selectMSigDBGeneSet(stag);
            catch ME            
                errordlg(ME.message);
                return;
            end
            ttxt = ctselected;
            if isempty(posg) || isempty(ctselected), return; end
            [y]=gui.e_cellscore(sce,posg);
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
            [~,cL,noanswer]=gui.i_reordergroups(thisc);
            if noanswer, return; end            
            colorit=true;
            cL=strrep(cL,'_','\_');
            thisc=strrep(thisc,'_','\_');
            gui.i_cascadeviolin(sce,Xt,thisc,glist, ...
                selitems{indx1},cL,colorit);
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

            %[glist]=gui.i_selectngenes(string(listitems));

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

    gui.i_violinplot(y,thisc,ttxt);

%     [~,cL]=grp2idx(thisc);
%     colorit=true;
%     f=figure('visible','off');
%     tb=uitoolbar(f);
%     pkg.i_addbutton2fig(tb,'off',{@i_savedata,y,thisc}, ...
%         'export.gif','Export data...');
%     pkg.i_addbutton2fig(tb,'off',{@i_testdata,y,thisc}, ...
%         'exportx.gif','ANOVA/T-test...');    
%     pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3}, ...
%         "powerpoint.gif",'Save Figure to PowerPoint File...');
%     pkg.i_addbutton2fig(tb,'off',@i_invertcolor, ...
%         "xpowerpoint.gif",'Switch BW/Color');
%     pkg.i_addbutton2fig(tb,'off',@i_reordersamples, ...
%         "xpowerpoint.gif",'Reorder Samples');   
% 
%     cL=strrep(cL,'_','\_');
%     thisc=strrep(thisc,'_','\_');
%     pkg.i_violinplot(y,thisc,colorit,cL);
%     title(strrep(ttxt,'_','\_'));
%     ylabel(selitems{indx1});
%     movegui(f,'center');
%     set(f,'visible','on');
%     
% %catch ME
% %    errordlg(ME.message);
% %end
% 
%     function i_invertcolor(~,~)
%         colorit=~colorit;
%         cla;
%         pkg.i_violinplot(y,thisc,colorit,cL);
%     end
% 
% 
%     function i_reordersamples(~,~)
%         [~,cL,noanswer]=gui.i_reordergroups(thisc);
%         if noanswer, return; end
%         cla
%         pkg.i_violinplot(y,thisc,colorit,cL);
%     end
% 
% end
% 
% function i_savedata(~,~,a,b)
%     T=table(a(:),b(:));    
%     T.Properties.VariableNames={'ScoreLevel','GroupID'};
%     T=sortrows(T,'ExprLevel','descend');
%     T=sortrows(T,'GroupID');
%     gui.i_exporttable(T,true);
% end
% 
% function i_testdata(~,~,y,grp)
%     if size(y,2)~=length(grp)
%         y=y.';
%     end
%     tbl=pkg.e_grptest(y,grp);
%     gui.i_exporttable(tbl,true);
% end
% 
