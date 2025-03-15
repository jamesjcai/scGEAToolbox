function [feedbk, extractedText] = e_GEOPageSummary(acc)
   
    feedbk = '';
    url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', acc);
    rawHtml = webread(url);

    %{
    if license('test', 'Text_Analytics_Toolbox')
        tree = htmlTree(htmlContent);
        extractedText = extractHTMLText(tree);
    else
        textOnly = regexprep(htmlContent, '<[^>]*>', ' ');
        textOnly = regexprep(textOnly, '\s+', ' ');
        extractedText = strtrim(textOnly);
    end
    %}

        % Remove script and style tags with their content
        rawHtml = regexprep(rawHtml, '<script.*?</script>', '', 'ignorecase');
        rawHtml = regexprep(rawHtml, '<style.*?</style>', '', 'ignorecase');
        
        % Remove HTML tags
        htmlContent = regexprep(rawHtml, '<[^>]*>', ' ');
        

% Remove HTML tags using regular expressions
text_content = regexprep(htmlContent, '<.*?>', '');

% Remove extra whitespace (multiple spaces, tabs, newlines)
text_content = regexprep(text_content, '\s+', ' ');

% Decode HTML entities like &nbsp; &amp; etc.
text_content = regexprep(text_content, '&nbsp;', ' ');
text_content = regexprep(text_content, '&amp;', '&');
text_content = regexprep(text_content, '&lt;', '<');
text_content = regexprep(text_content, '&gt;', '>');
text_content = regexprep(text_content, '&quot;', '"');

% Trim leading/trailing whitespace
extractedText = strtrim(text_content);




    if isempty(which('ollamaChat'))
        error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
    end
    % https://www.mathworks.com/matlabcentral/fileexchange/163796-large-language-models-llms-with-matlab/
    % Add-On “Large Language Models (LLMs) with MATLAB”.
    
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
    
    try
        chat = ollamaChat(providermodel{2}, TimeOut = 1200);
        prompt1 = "Please summerize the following text downloaded from the GEO database. ";
        prompt2 = "Here is the text: " + extractedText;
        feedbk = generate(chat, prompt1 + prompt2);
        feedbk = regexprep(feedbk, '<think>.*?</think>', '');
    catch ME
        warning(ME.message);
    end

    %{
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
    close(doc);
    rptview(outfile, 'docx');
    %}
end