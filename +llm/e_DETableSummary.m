function [done, outfile] = e_DETableSummary(TbpUp, TmfUp, ...
                            TbpDn, TmfDn, infotagstr)
    
    done = false;
    outfile = [];
    
    % Steps to test this function:
    %
    % (1) Download and install Ollama https://ollama.com/
    % (2) Use Ollama to pull a model, e.g., Ollama pull deepseek-r1
    % (3) In MATLAB, install the Add-On, Large Language Models (LLMs) with MATLAB
    % (4) In SCGEATOOL, click menu Options -> Set LLM Provider & Model...
    % (5) In SCGEATOOL, click menu Analysis -> All (DE, DV, DP) Analysis in Cell Type Batch Mode...
    % (6) In MATLAB Command Window, run llm.function
    
    if isempty(which('ollamaChat'))
        error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
    end
    % https://www.mathworks.com/matlabcentral/fileexchange/163796-large-language-models-llms-with-matlab/
    % Add-On “Large Language Models (LLMs) with MATLAB”.
    
    if nargin<4, infotagstr = pkg.i_randinfostr; end
        
    wrkdir = getpref('scgeatoolbox', 'externalwrkpath');
    if isempty(wrkdir), return; end
        
    preftagname = 'llmodelprovider';
    s = getpref('scgeatoolbox', preftagname);
    providermodel = strsplit(s,':');
    assert(strcmp(providermodel{1}, 'Ollama'))
    
    try
        aa=webread("http://localhost:11434");
        disp(aa)
    catch ME
        disp('Ollama is not running.');
        return;
    end
    
    warning off
    s_up = ''; 
    s_dn = '';
    if ~isempty(TbpUp) || ~isempty(TmfUp)
        T1 = []; T2 = [];
        if ~isempty(TbpUp)
            T1 = TbpUp(:, [3 7]);
        end
        if ~isempty(TmfUp)
            T2 = TmfUp(:, [3 7]);
        end        
        T = [T1; T2];
        T.Properties.VariableNames={'Gene Ontology Term','Genes'};
        
        % T.Genes = strrep(T.Genes,',', ', ');
        C = T.Genes;
        C_new = cellfun(@(x) strrep(x, ',', ', '), C, 'UniformOutput', false);
        T.Genes = C_new;
    
        s_up = jsonencode(table2struct(T));
    end
    
    if ~isempty(TbpDn) || ~isempty(TmfDn)
        T1 = []; T2 = [];
        if ~isempty(TbpDn)
            T1 = TbpDn(:, [3 7]);
        end
        if ~isempty(TmfDn)
            T2 = TmfDn(:, [3 7]);
        end
        T = [T1; T2];
        T.Properties.VariableNames={'Gene Ontology Term','Genes'};
        C = T.Genes;
        C_new = cellfun(@(x) strrep(x, ',', ', '), C, 'UniformOutput', false);
        T.Genes = C_new;
        s_dn = jsonencode(table2struct(T));
    end
    warning on
    
    if isempty(s_up) && isempty(s_dn), return; end
       
    chat = ollamaChat(providermodel{2}, TimeOut = 1200);
    prompt1 = "You are a researcher working on gene function enrichment analysis try to find biological processes and molecular functions of a given gene sets. You are using Enrichr to do this. " + ...
        "Enrichr is a gene function enrichment analysis tool. You will be given the output of Enrichr analysis below, which contains a list of gene ontology (GO) terms and associated genes. The GO terms are a mix of enriched terms of biological processes and molecular functions. Please summarize the results in text. " + ...
        "Please provide an analysis of the output, highlighting key biological processes and molecular functions, along with their associated genes. " + ...
        "Please write an executive summary to report the results of your analysis. ";
    prompt2 = "Here is the output of Enrichr: " + s_up;
    feedbk_up = generate(chat, prompt1 + prompt2);
    
    prompt2 = "Here is the output of Enrichr: " + s_dn;
    feedbk_dn = generate(chat, prompt1 + prompt2);
    
    import mlreportgen.dom.*
    
    outfile = fullfile(wrkdir, "Res_"+matlab.lang.makeValidName(infotagstr));
    doc = Document(outfile, 'docx');
    open(doc);
    
    para = Paragraph(sprintf('%s Up-regulation', infotagstr));
    para.Style = {Bold(true), FontSize('14pt'), Color('blue')};
    append(doc, para);
    
    feedbk_up = regexprep(feedbk_up, '<think>.*?</think>', '');
    para = Paragraph(feedbk_up);
    append(doc, para);
    
    para = Paragraph(sprintf('%s Down-regulation', infotagstr));
    para.Style = {Bold(true), FontSize('14pt'), Color('blue')};
    append(doc, para);
    
    feedbk_dn = regexprep(feedbk_dn, '<think>.*?</think>', '');
    para = Paragraph(feedbk_dn);
    append(doc, para);
    
    close(doc);
    %rptview(outfile, 'docx');
    done = true;

end