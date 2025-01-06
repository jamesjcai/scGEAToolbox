function callback_RunEnrichr(src, ~, predefinedlist, enrichrtype, ...
    backgroundlist, outfiletag)

if nargin < 6, outfiletag = ""; end
if nargin < 5
    askbackground = true;
    backgroundlist = []; 
else
    askbackground = false;
end
if nargin < 4, enrichrtype = []; end
if nargin < 3, predefinedlist = []; end

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    gsorted = natsort(sce.g);

    rng("shuffle");
    n = length(gsorted);
    if isempty(predefinedlist)
        error('No input gene list.');
        % ingenelist = gui.i_inputgenelist(gsorted(randperm(n, min([200, length(gsorted)]))));
    else
        ingenelist = gui.i_inputgenelist(predefinedlist);
    end
    if isempty(ingenelist), return; end

    if askbackground
        answer = questdlg('Add background list?','');
        switch answer
            case 'Yes'
                [idx] = gui.i_selmultidlg(sce.g, sce.g, FigureHandle);
                if isempty(idx), return; end
                if idx == 0, return; end
                backgroundlist = sce.g(idx);
            case 'No'
                backgroundlist = [];
            case 'Cancel'
                return;
        end
    end

    if isempty(enrichrtype)
        answer1 = questdlg("Select the type of Enrichr application.","", ...
             "Web-based", "API-based", "API-based");
        enrichrtype = answer1;
    end

switch enrichrtype 
    case "API-based"
        % do nothing here
    case "Web-based"
        fw = gui.gui_waitbar([], false, 'Sending genes to web browser...');
        % gui.i_enrichtest(genelist, backgroundlist, numel(genelist));
            if ~isempty(backgroundlist)
                run.web_Enrichr_bkg(ingenelist, backgroundlist, numel(ingenelist));
            else
                run.web_Enrichr(ingenelist, numel(ingenelist));
            end
        gui.gui_waitbar(fw, false, 'Check web browser & submit genes to Enrichr.');
        return;   
    otherwise
        return;
end



    [genesets] = in_selDataSources;
    if isempty(genesets), return; end

    definput = {'5', '0.1'};
    prompt = {'Min # of overlapping genes:', ...
              'P-value cutoff:'};
    dlgtitle = 'Enrichr Result Filter';
    dims = [1, 80];
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer), return; end
    try
        minugenes = str2double(answer{1});
        pvaluecut = str2double(answer{2});
        assert((minugenes > 0) && (minugenes < 100));
        assert((pvaluecut >= 0.0) && (pvaluecut <= 1.0));
    catch
        errordlg('Invalid input.');
        return;
    end
    

    fw = gui.gui_waitbar;
    Tlist = run.ml_Enrichr(ingenelist, backgroundlist, genesets,...
                           minugenes, pvaluecut);

    T=table;
    for k = 1:height(Tlist)
        if ~isempty(Tlist{k}) && istable(Tlist{k})
            T = [T; Tlist{k}];
        end
    end

    gui.gui_waitbar(fw);
    
    [~, ~] = gui.i_exporttable(T, true, 'Tenrichrres', ...
        sprintf('Enrichr_Results_%s', outfiletag));

    % gui.i_viewtable(T, FigureHandle);


    function [genesets] = in_selDataSources
        genesets = [];
        dsv = pkg.i_get_enrichr_libraries;
        % if ~ispref('scgeatoolbox', 'enrichrlibraries')
        %     setpref('scgeatoolbox', 'enrichrlibraries', ["GO_Biological_Process_2023", ...
        %                 "GO_Molecular_Function_2023"]);
        % end
        enrichrlibraries = getpref('scgeatoolbox', 'enrichrlibraries', ...
                                ["GO_Biological_Process_2023", ...
                                    "GO_Molecular_Function_2023"]);

        [idx1] = gui.i_selmultidlg(dsv, enrichrlibraries, FigureHandle);
        if isempty(idx1), return; end
        if idx1 == 0, return; end
        genesets = dsv(idx1);
        setpref('scgeatoolbox', 'enrichrlibraries', genesets);    
    end

end
