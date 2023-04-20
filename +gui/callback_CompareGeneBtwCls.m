function callback_CompareGeneBtwCls(src,~)
    answer = questdlg('This function generates violinplot to show differences between cell groups. Continue?','');
    if ~strcmp(answer,'Yes'), return; end    

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc]=gui.i_select1class(sce);
    if isempty(thisc), return; end
    

    if length(unique(thisc))==1
        answer=questdlg("All cells are in the same group. No comparison will be made. Continue?", ...
            "",'Yes','No','Cancel','No');
        switch answer
            case 'Yes'
            otherwise
                return;
        end
    else
        [ci,cLi]=grp2idx(thisc);
        listitems=natsort(string(cLi));
        n=length(listitems);
        [indxx,tfx] = listdlg('PromptString',{'Select two groups:'},...
            'SelectionMode','multiple',...
            'ListString',listitems,...
            'InitialValue',1:n);
        if tfx==1
            [y1,idx1]=ismember(listitems(indxx),cLi);
            assert(all(y1));
            idx2=ismember(ci,idx1);
            sce=sce.selectcells(idx2);
            thisc=thisc(idx2);
        else
            return;
        end    

    end

% selitems={'Expression of Gene', ...
%     'TF Activity Score [PMID:33135076]',...
%     'TF Targets Expression Score',...
%     'Differentiation Potency [PMID:33244588]',...
%     'MSigDB Signature Score',...
%     '--------------------------------',...
%     'Predefined Cell Score',...
%     'Define New Score...',...
%     '--------------------------------',...
%     'Library Size','Other Attribute'};
% [indx1,tf1]=listdlg('PromptString',...
%     'Select a metric for comparison.',...
%     'SelectionMode','single','ListString',selitems, ...
%     'ListSize',[200,300]);
% if tf1~=1, return; end
% selecteditem=selitems{indx1};

[selecteditem] = gui.i_selectcellscore;
%try
    switch selecteditem
        %case 'Global Coordination Level (GCL) [PMID:33139959]'

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
        case 'Differentiation Potency [PMID:33244588]'
                answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                [yes,speciesid]=ismember(lower(answer2),{'human','mouse'});
                if ~yes, return; end
                y=sc_potency(sce.X,sce.g,speciesid);
                ttxt='Differentiation Potency';

            % [a]=contains(sce.list_cell_attributes(1:2:end),'cell_potency');
            % if ~any(a)
            %     answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
            %     [yes,specisid]=ismember(lower(answer2),{'human','mouse'});
            %     if ~yes, return; end
            %     sce=sce.estimatepotency(specisid);
            % end
            % [yes,idx]=ismember({'cell_potency'},sce.list_cell_attributes(1:2:end));
            % if yes
            %     y=sce.list_cell_attributes{idx*2};                
            %     ttxt='Differentiation Potency';
            % else                
            %     return;
            % end

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
            if isstring(thisc)
                thisc=strrep(thisc,'_','\_');
            end
            gui.i_cascadeviolin(sce,Xt,thisc,glist, ...
                'Expression Level',cL,colorit);
            return;
        case 'Predefined Cell Score'
            [~,T]=pkg.e_cellscores(sce.X,sce.g,0);
            listitems=T.ScoreType;
            [indx2,tf2] = listdlg('PromptString','Select Class',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end
            %fw=gui.gui_waitbar;
            y=pkg.e_cellscores(sce.X,sce.g,indx2);
            ttxt=T.ScoreType(indx2);
            %gui.gui_waitbar(fw);
        case 'Other Attribute...'
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
        case {'TF Activity Score [PMID:33135076] 🐢',...
                'TF Targets Expression Score'}
            [~,T]=pkg.e_tfactivityscores(sce.X,sce.g,0);
            listitems=unique(T.tf);

            %[glist]=gui.i_selectngenes(string(listitems));

            [indx2,tf2] = listdlg('PromptString','Select Class',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end
            species=gui.i_selectspecies(2);
            if isempty(species), return; end
            
%[cs]=gui.e_cellscore(sce,posg);
%     answer = questdlg('Select algorithm:',...
%     'Select Method', ...
%     'UCell [PMID:34285779]','AddModuleScore/Seurat', ...
%     'UCell [PMID:34285779]');
%     switch answer
%         case 'AddModuleScore/Seurat'   

            if strcmp(selecteditem,'TF Targets Expression Score')
                methodid=1;
            else
                methodid=4;
            end
            if methodid~=4
            fw=gui.gui_waitbar;
            end
                [cs,tflist]=sc_tfactivity(sce.X,sce.g,[],species,methodid);
                idx=find(tflist==string(listitems{indx2}));
                assert(length(idx)==1)
                [y]=cs(idx,:);
                ttxt=listitems{indx2};
            if methodid~=4
            gui.gui_waitbar(fw);
            end

        case 'TF Targets Expression Score 2'
                species=gui.i_selectspecies(2);
                if isempty(species), return; end            
                [posg,ctselected]=gui.i_selectTFTargetSet(species);
                [y]=gui.e_cellscore(sce,posg);
                %ttxt=listitems{indx2};
                ttxt = strcat(ctselected, ' activity');

        otherwise
            return;
    end
% assignin('base','y',y);
% assignin('base','thisc',thisc);

    gui.i_violinplot(y,thisc,ttxt);

end