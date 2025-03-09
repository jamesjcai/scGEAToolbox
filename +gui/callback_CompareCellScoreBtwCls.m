function [needupdatesce] = callback_CompareCellScoreBtwCls(src, ~)

%{
1. Input:
   - `sce`: An object containing the cellular data.
   - `glist`: A list of gene/region names (optional).
   - `thisc`: The current cell group being analyzed.
   - `selecteditem`: Determines which type of score to calculate or display.

2. Output:
    - Displays a graph showing the calculated scores along with labels for each dataset.
%}

[FigureHandle, sce, isui] = gui.gui_getfigsce(src);
needupdatesce = false;

aa = 'Yes, compare scores (violinplot)';
bb = 'No, just show values (heatmap)';
    answer2 = questdlg(sprintf(['This function will calculates a score for each cell. After the scores are calculated, do you want to ', ...
                'compare score values between different cell groups?']), '', ...
                bb, aa, bb);
            switch answer2
                case aa
                    showcomparision = true;
                case bb
                    showcomparision = false;
                otherwise
                    return;
            end

            if ~showcomparision
                thisc = ones(sce.NumCells, 1);
            else
                allowunique = false;
                [thisc] = gui.i_select1class(sce, allowunique);
                if isempty(thisc), return; end
                if isscalar(unique(thisc))
                    answer = questdlg("All cells are in the same group. No comparison will be made. Continue?", ...
                        "", 'Yes', 'No', 'Cancel', 'No');
                    switch answer
                        case 'Yes'
                        otherwise
                            return;
                    end
                else    % length(unique(thisc)) ~= 1
                    [ci, cLi] = grp2idx(thisc);
                    listitems = natsort(string(cLi));
                    n = length(listitems);
                    [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
                        'SelectionMode', 'multiple', ...
                        'ListString', listitems, ...
                        'InitialValue', 1:n, ...
                        'ListSize', [220 300]);
                    if tfx == 1
                        [y1, idx1] = ismember(listitems(indxx), cLi);
                        assert(all(y1));
                        idx2 = ismember(ci, idx1);
                        sce = sce.selectcells(idx2);
                        thisc = thisc(idx2);
                    else
                        return;
                    end
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
            %     'ListSize',[220, 300]);
            % if tf1~=1, return; end
            % selecteditem=selitems{indx1};

            % ------------------------------------------------

            [selecteditem, speciestag] = gui.i_selgenesetcollection;
            if isempty(selecteditem), return; end
            %try

            switch selecteditem
                    %case 'Global Coordination Level (GCL) [PMID:33139959]'

                case 'Define a New Score...'
                    ttxt = 'Customized Score';
                    newcstype = inputdlg('Name of the new cell score:', '', ...
                            [1, 50], {ttxt});
                    if isempty(newcstype), return; end
                    ttxt = newcstype;
                    [posg] = gui.i_selectngenes(sce.g, [], FigureHandle);
                    if isempty(posg)
                        helpdlg('No feature genes selected.', '')
                        return;
                    end
                    [y] = gui.e_cellscore(sce, posg);
                    answer = questdlg('Save score to cell attribute list?', ...
                        '','Yes, save','No, skip','Cancel','Yes, save');
                    switch answer
                        case 'Yes, save'
                            ttxt = string(ttxt);
                            if any(strcmp(ttxt, sce.list_cell_attributes(1:2:end)))                            
                                fprintf("Duplicate found: %s is renamed to ", ttxt); 
                                a = matlab.lang.makeUniqueStrings([string(sce.list_cell_attributes(1:2:end)) ttxt]);
                                ttxt = a(end);
                                fprintf("%s.\n", ttxt);
                            else
                                sce.list_cell_attributes = [sce.list_cell_attributes, ...
                                    {ttxt, y(:)}];
                            end
                        case 'No, skip'
                        otherwise
                            return;
                    end
                case 'MSigDB Molecular Signatures'
                    % speciestag = gui.i_selectspecies(2, true);
                    % if isempty(speciestag), return; end
                    try
                        [posg, ctselected] = gui.i_selectMSigDBGeneSets(speciestag);
                    catch ME
                        errordlg(ME.message);
                        return;
                    end
                    ttxt = ctselected;
                    if isempty(posg) || isempty(ctselected), return; end

                    n = length(posg);
                    y = cell(n,1);
                    [~, methodid] = gui.i_pickscoremethod;
                    if isempty(methodid), return; end

                    fw = gui.gui_waitbar;
                    for k = 1:n
                        y{k} = gui.e_cellscore(sce, posg{k}, ...
                            methodid, false);
                    end
                    gui.gui_waitbar(fw);
                    

                case 'PanglaoDB Cell Type Markers'
                    if isempty(speciestag)
                        speciestag = gui.i_selectspecies(2, true);
                    end
                    if isempty(speciestag), return; end

                    oldpth = pwd;
                    pw1 = fileparts(mfilename('fullpath'));
                    pth = fullfile(pw1, '..', '+run', 'thirdparty', 'alona_panglaodb');
                    cd(pth);
                    switch lower(speciestag)
                        case {'human', 'hs'}
                            speciestag_short = 'hs';
                        case {'mouse', 'mm'}
                            speciestag_short = 'mm';
                        otherwise
                            speciestag_short = 'hs';
                    end


                    markerfile = sprintf('marker_%s.mat', speciestag_short);
                    if exist(markerfile, 'file')
                        load(markerfile, 'Tm');
                    else
                        % Tw=readtable(sprintf('markerweight_%s.txt',stag));
                        Tm = readtable(sprintf('markerlist_%s.txt', speciestag_short), ...
                            'ReadVariableNames', false, 'Delimiter', '\t');
                        % save(markerfile,'Tw','Tm');
                    end
                    cd(oldpth);

                    ctlist = string(Tm.Var1);
                    listitems = sort(ctlist);
                    [indx, tf] = listdlg('PromptString', ...
                        {'Select Class'}, ...
                        'SelectionMode', 'single', 'ListString', listitems, ...
                        'ListSize', [220, 300]);
                    if ~tf == 1, return; end
                    ctselected = listitems(indx);
                    % idx=find(matches(ctlist,ctselected));
                    idx = matches(ctlist, ctselected);
                    ctmarkers = Tm.Var2{idx};
                    posg = string(strsplit(ctmarkers, ','));
                    posg(strlength(posg) == 0) = [];
                    ttxt = ctselected;
                    if isempty(posg) || isempty(ctselected), return; end
                    
                    [y] = gui.e_cellscore(sce, posg);
                    

                % case 'Differentiation Potency [PMID:33244588]'
                % 
                %     if ~gui.gui_showrefinfo('Differentiation Potency [PMID:33244588]'), return; end
                % 
                %     % answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                %     % [yes,speciesid]=ismember(lower(answer2),{'human','mouse'});
                %     % if ~yes, return; end
                %     speciestag = gui.i_selectspecies(2);
                % 
                %     fw = gui.gui_waitbar;
                %     y = sc_potency(sce.X, sce.g, speciestag);
                %     if sum(strcmp('cell_potency', sce.list_cell_attributes)) == 0
                %         sce.list_cell_attributes = [sce.list_cell_attributes, ...
                %             {'cell_potency', y(:)}];                  
                %     end
                %     if ~showcomparision
                %         needupdatesce = true;
                %         guidata(FigureHandle, sce);
                %     end
                % 
                %     ttxt = 'Differentiation Potency';
                %     posg = [];
                %     gui.gui_waitbar(fw);

                % [a]=contains(sce.list_cell_attributes(1:2:end),'cell_potency');
                % if ~any(a)
                %     answer2=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                %     [yes,specisid]=ismember(lower(answer2),{'human','mouse'});
                %     if ~yes, return; end
                %     sce=sce.estimatepotency(specisid);
                % end
                % [yes, idx]=ismember({'cell_potency'},sce.list_cell_attributes(1:2:end));
                % if yes
                %     y=sce.list_cell_attributes{idx*2};
                %     ttxt='Differentiation Potency';
                % else
                %     return;
                % end

                % case 'Library Size of Cells'
                %     y = sum(sce.X);
                %     ttxt = 'Library Size';
                %     posg = [];
                %     if sum(strcmp('library_size', sce.list_cell_attributes)) == 0
                %         sce.list_cell_attributes = [sce.list_cell_attributes, ...
                %             {'library_size', y(:)}];                  
                %     end                    
                % 
                % case 'Expression of Individual Genes'
                %     [glist] = gui.i_selectngenes(sce,[],FigureHandle);
                %     if isempty(glist)
                %         helpdlg('No gene selected.', '');
                %         return;
                %     end
                %     [Xt] = gui.i_transformx(sce.X);
                %     [~, cL, noanswer] = gui.i_reordergroups(thisc);
                %     if noanswer, return; end
                %     colorit = true;
                %     cL = strrep(cL, '_', '\_');
                %     if isstring(thisc)
                %         thisc = strrep(thisc, '_', '\_');
                %     end
                %     gui.i_cascadeviolin(sce, Xt, thisc, glist, ...
                %         'Expression Level', cL, colorit);
                %     return;

                case 'Predefined Custom Gene Sets'
                    if ~gui.gui_showrefinfo('Predefined Cell Score'), return; end
                    [~, T] = pkg.e_cellscores(sce.X, sce.g, 0);
                    listitems = T.ScoreType;
                    [indx2, tf2] = listdlg('PromptString', 'Select Score', ...
                        'SelectionMode', 'multiple', 'ListString', ...
                        listitems, 'ListSize', [260, 300]);
                    if tf2 ~= 1, return; end

                    n = length(indx2);
                    y = cell(n,1); ttxt=cell(n,1);
                    [~, methodid] = gui.i_pickscoremethod;

                    fw=gui.gui_waitbar;
                    for k = 1:n
                        [y{k}, ~, posg] = pkg.e_cellscores(sce.X, sce.g, ...
                            indx2(k), methodid, false);
                        ttxt{k} = T.ScoreType(indx2(k));
                    end
                    gui.gui_waitbar(fw);
                    % case 'Other Cell Attribute...'
                    %     [y, clabel, ~, newpickclabel] = gui.i_select1state(sce, true);
                    %     if isempty(y)
                    %         helpdlg('No cell attribute is available.');
                    %         return;
                    %     end
                    %     if ~isempty(newpickclabel)
                    %         ttxt = newpickclabel;
                    %     else
                    %         ttxt = clabel;
                    %     end
                    %     posg = [];
                case {'TF Activity Score [PMID:33135076] 🐢', ...
                      'DoRothEA TF Targets'}
                    [T] = pkg.e_gettflist;
                    listitems = unique(T.tf);
                   
                    [indx2, tf2] = listdlg('PromptString', 'Select a transcription factor (TF)', ...
                        'SelectionMode', 'single', 'ListString', ...
                        listitems, 'ListSize', [220, 300]);
                    if tf2 ~= 1, return; end
                    %species = gui.i_selectspecies(2);
                    %if isempty(species), return; end

                    methodid = 4;

                    % if strcmp(selecteditem, 'DoRothEA TF Targets')
                    %     methodid = 4;  % UCell method
                    %     % disp('Using the UCell method...');
                    % else
                    %     methodid = 4;
                    % end
                    
                    if methodid ~= 4, fw = gui.gui_waitbar; end
                    try
                        [cs, tflist] = sc_tfactivity(sce.X, sce.g, [], ...
                            speciestag, methodid);
                    catch ME
                        if methodid ~= 4, gui.gui_waitbar(fw,true); end
                        errordlg(ME.message);
                        return;
                    end
                    idx = find(upper(tflist) == upper(string(listitems{indx2})));
                    assert(isscalar(idx))

                    [y] = cs(idx, :);
                    ttxt = listitems{indx2};
                    posg = [];
                    if methodid ~= 4, gui.gui_waitbar(fw); end
                    %         case 'TF Targets Expression Score 2'
                    %                 species=gui.i_selectspecies(2);
                    %                 if isempty(species), return; end
                    %                 [posg,ctselected]=gui.i_selectTFTargetSet(species);
                    %                 [y]=gui.e_cellscore(sce,posg);
                    %                 %ttxt=listitems{indx2};
                    %                 ttxt = strcat(ctselected, ' activity');

                otherwise
                    return;
            end
            %catch ME
            %    %errordlg(ME.message);
            %    return;
            %end

            if isempty(y), return; end

            % assignin('base','y',y);
            % assignin('base','thisc',thisc);

            if ~exist("posg", "var"), posg = []; end
            if ~exist("ttxt", "var"), ttxt = []; end

            if showcomparision
                %if iscell(y)
                    gui.sc_uitabgrpfig_vioplot(y, ttxt, thisc, FigureHandle);
                %else
%                    gui.i_violinplot(y, thisc, ttxt, true, [], posg, FigureHandle);
%                    xlabel('Cell group');
%                    ylabel('Cellular score');
%                end
            else
                %     [methodid]=gui.i_pickscatterstem('Scatter+Stem');
                %     if isempty(methodid), return; end
                %         f=gui.i_cascadefig(sce,glist(k),axx,bxx,k,methodid);
                %     [h1]=sc_scattermarker(sce.X,sce.g,...
                %                  sce.s,g,methodid);
                if iscell(y)
                    % t=array2table(cell2mat({rand(10,1),rand(10,1),rand(10,1)}),'VariableNames',{'aa','bb','cc'});                    
                    % assignin("base",'y',y);
                    % assignin("base",'ttxt',ttxt);
                    % assignin("base",'k',k);                  

                    if length(y)>1
                        pause(1);
                        gui.i_scoreheatmap(cell2mat(y').', ttxt, sce, FigureHandle);
                    elseif isscalar(y)
                        gui.i_stemscatterfig(sce, y{1}, posg, ttxt{1}, FigureHandle);
                    end
                        % gui.sc_uitabgrpfig_expplot(y, markerlist, sce.s, FigureHandle, [axx, bxx]);
                else

                    if sum(strcmp(ttxt, sce.list_cell_attributes)) == 0
                        sce.list_cell_attributes = [sce.list_cell_attributes, ...
                            {ttxt, y(:)}];
                    end
                    needupdatesce = true;
                    guidata(FigureHandle, sce);
                    gui.i_stemscatterfig(sce, y, posg, ttxt, FigureHandle);
                end
            end
            % guidata(FigureHandle, sce);

end
