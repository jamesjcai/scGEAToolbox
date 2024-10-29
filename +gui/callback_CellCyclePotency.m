function [needupdate] = callback_CellCyclePotency(src, ~, typeid)

    if nargin < 3, typeid = 1; end

    needupdate = false;
    needestimt = false;

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    
    switch typeid
        case 1
            answer = questdlg('This function assigns cell cycle phase to each cell, continue?', '');
        case 2
            answer = questdlg('This function assigns differentiation potency [PMID:33244588] to each cell, continue?', '');
        case 3
            answer = questdlg('This function calculates stemness index [PMID:29625051] for each cell, continue?', '');
        case 4
            answer = questdlg('This function calculates the expression ratio of dissociation-associated genes [PMID:34020534] for each cell, continue?', '');
        case 5
            answer = questdlg('This function predicts tumor (aneuploid) and normal (diploid) cells using copykat [PMID:33462507], continue?', '');
    end
    if ~strcmp(answer, 'Yes'), return; end


    switch typeid
        case 1
            if isempty(sce.c_cell_cycle_tx) || all(strcmp(unique(sce.c_cell_cycle_tx), "undetermined"))
                needestimt = true;
            else
                answer1 = questdlg('Use existing cell cycle estimation or re-compute new estimation?', ...
                    '', 'Use existing', 'Re-compute', 'Cancel', 'Use existing');
                switch answer1
                    case 'Re-compute'
                        needestimt = true;
                    case 'Cancel'
                        return;
                end
            end
            if needestimt
                fw = gui.gui_waitbar;
                sce = sce.estimatecellcycle(true, 1);
                needupdate = true;
                gui.gui_waitbar(fw);
                guidata(FigureHandle, sce);
                uiwait(helpdlg('Cell cycle phase (c_cell_cycle_tx) added.', ''));
            end
            y = sce.c_cell_cycle_tx;
            attribtag = "cell_cycle";
            uiwait(helpdlg('To see the result, use View -> Cell State (Ctrl + T). Then select "Cell Cycle Phase"',''));
            
        case 2
            attribtag = "cell_potency";
            y = in_aaa(attribtag);
        case 3
            attribtag = "stemness_index";
            y =in_aaa(attribtag);
        case 4
            attribtag = 'dissocation_ratio';
            y = in_aaa(attribtag);
        case 5
            attribtag = 'copykat_prediction';
            y = in_aaa(attribtag);
    end


    function [s] = in_aaa(attribtag)        
        if ~ismember(attribtag, sce.list_cell_attributes(1:2:end))
            needestimt = true;
        else
            answer1 = questdlg(sprintf('Use existing %s estimation or re-compute new estimation?', ...
                attribtag), '', 'Use existing', 'Re-compute', 'Cancel', 'Use existing');
            switch answer1
                case 'Re-compute'
                    needestimt = true;
                case 'Cancel'
                    return;
            end
        end
        if needestimt
            switch attribtag
                case 'cell_potency'
                    speciestag = gui.i_selectspecies(2, false);
                    if isempty(speciestag), return; end
            end

            fw = gui.gui_waitbar;            
            switch attribtag
                case 'cell_potency'
                    sce = sce.estimatepotency(speciestag);
                    needupdate = true;
                    [yesx, idx] = ismember('cell_potency', sce.list_cell_attributes(1:2:end));
                    assert(yesx);
                    s =  sce.list_cell_attributes{idx+1};
                case 'stemness_index'
                    s = sc_stemness(sce.X, sce.g);
                case 'dissocation_ratio'
                    s = pkg.sc_dissratio(sce.X, sce.g, true);
                case 'copykat_prediction'
                    s = run.r_copykat(sce);
                otherwise
                    error('Invalid attribtag');
            end
            [yesx, idx] = ismember(attribtag, sce.list_cell_attributes(1:2:end));
            if yesx
                sce.list_cell_attributes{idx*2} = s;
            else
                sce.list_cell_attributes = [sce.list_cell_attributes, ...
                    {attribtag, s}];
            end
            gui.gui_waitbar(fw);
            guidata(FigureHandle, sce);
            needupdate = true;
            uiwait(helpdlg(sprintf('%s added.', attribtag), ''));
        end
        uiwait(helpdlg(sprintf('To see the result, use View -> Cell State (Ctrl + T). Then select "%s"', ...
                attribtag),''));
    end

end