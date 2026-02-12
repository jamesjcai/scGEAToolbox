function callback_RunGeneAgent(src, ~, predefinedlist)

% [PMID:40721871]

if nargin < 3, predefinedlist = []; end

if ~isempty(predefinedlist)
    [FigureHandle] = gui.gui_getfigsce(src);
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end


    if ~pkg.i_license
        gui.myErrordlg(FigureHandle, ...
            "This function requires passkey validation. You can " + ...
            "validate your passkey by selecting Help → Validate " + ...
            " Passkey from the menu.");
        return;
    end


    extprogname = 'geneagentwork';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, ...
        preftagname, FigureHandle);
    if isempty(wkdir), return; end
    olddir = pwd;
    if isfolder(wkdir), cd(wkdir); end
    
    if isempty(predefinedlist)
        gsorted = natsort(sce.g);
        rng("shuffle");
        n = length(gsorted);
        if isempty(gsorted)
            ingenelist = gui.i_inputgenelist(sprintf("ERBB2\nERBB4\nFGFR2\nFGFR4\nHRAS\nKRAS"));
        else
            ingenelist = gui.i_inputgenelist(gsorted(randperm(n, ...
                min([7, length(gsorted)]))));
        end
    else
        ingenelist = predefinedlist;
    end

    if isempty(ingenelist) || all(strlength(ingenelist) < 1), return; end
    fw = gui.myWaitbar(FigureHandle);
    ingenelist = sprintf("%s,", ingenelist);
    ingenelist = extractBefore(ingenelist, strlength(ingenelist));

    try
        [~, retrieveurl] = run.ml_geneagent(ingenelist);        
    catch ME
        cd(olddir);
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message);
        return;
    end
    gui.myWaitbar(FigureHandle, fw);

    if strcmp('Yes', gui.myQuestdlg(FigureHandle, "Wait until the analysis is complete, then generate the report."))
        
           gui.myHelpdlg(FigureHandle, ...
            'Wait until the web application is done, click OK to proceed.', ...
            'Process Status');
        
        % The code execution will pause here until user clicks Continue
        % if strcmp(selection, 'Continue')
            options = weboptions('Timeout', 30);
            out = webread(retrieveurl, options);
            % Process the output and generate the report
            if isstruct(out)
                in_generateAIReport(out);
            else
                warning('Output is not a struct. Report generation skipped.');
            end
        % end
    end
    cd(olddir);
end



function in_generateAIReport(s)
    str = formattedDisplayText(s);
    writelines(str, 'GeneAgent_Report.txt');
    % type 'struct_report.txt';    
    import mlreportgen.dom.*    
    % Assume your struct with 10 text paragraphs:
    
    fields = fieldnames(s);
    % for k = 1:numel(fields)
    %     s.(fields{k}) = sprintf('Long paragraph for %s...', fields{k}); 
    % end
    
    % Create a Word report:
    doc = Document('GeneAgent_Report','docx');
    
    % Optional: Add a title
    append(doc, Heading1('GeneAgent Report'));
    
    % Loop through struct fields
    for k = 1:numel(fields)
        fname = fields{k};
        txt = s.(fname);
        
        % Add a subheading for field name
        append(doc, Heading2(sprintf('%s', fname)));
        
        % Create formatted paragraph
        p = Paragraph(txt);
        p.Style = {FontFamily('Times New Roman'), FontSize('11pt'), ...
            OuterMargin('0.25in','0in','0in','12pt')};
        append(doc, p);
    end
    % Finalize and open the report
    close(doc);
    rptview(doc.OutputPath);
end