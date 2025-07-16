function geoInfo = extractGEOInfo(url_or_accession)
    % EXTRACTGEOINFO Extract and analyze information from GEO database records
    %
    % SYNTAX:
    %   geoInfo = extractGEOInfo(url_or_accession)
    %
    % INPUT:
    %   url_or_accession - Either a full GEO URL or just the accession number (e.g., 'GSM7855468')
    %
    % OUTPUT:
    %   geoInfo - Structure containing extracted and analyzed information
    %
    % EXAMPLE:
    %   info = extractGEOInfo('GSM7855468');
    %   info = extractGEOInfo('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7855468');
    
    % Ensure Text Analytics Toolbox is available
    if ~license('test', 'Text_Analytics_Toolbox')
        error('Text Analytics Toolbox is required for this function');
    end
    
    % Parse input to get URL
    if startsWith(url_or_accession, 'http')
        url = url_or_accession;
        % Extract accession from URL
        tokens = regexp(url, 'acc=([^&]+)', 'tokens');
        if ~isempty(tokens)
            accession = tokens{1}{1};
        else
            accession = 'Unknown';
        end
    else
        accession = url_or_accession;
        url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s', accession);
    end
    
    try
        % Fetch the web page content
        fprintf('Fetching GEO record for %s...\n', accession);
        htmlContent = webread(url);
        
        % Extract plain text from HTML
        plainText = llm.extractplaintext(htmlContent);
        
        % Initialize output structure
        geoInfo = struct();
        geoInfo.accession = accession;
        geoInfo.url = url;
        geoInfo.rawText = plainText;
        
        % Extract structured information using regex patterns
        geoInfo.structured = extractStructuredInfo(plainText);
        
        % Use LLM to analyze and extract key information
        geoInfo.llm_analysis = analyzWithLLM(plainText, accession);
        
        % Use Text Analytics Toolbox for additional analysis
        geoInfo.text_analysis = performTextAnalysis(plainText);
        
        % Generate summary
        geoInfo.summary = generateSummary(geoInfo);
        
        fprintf('Successfully extracted information for %s\n', accession);
        
    catch ME
        warning(ME.identifier, 'Error processing GEO record: %s', ME.message);
        geoInfo = struct();
        geoInfo.accession = accession;
        geoInfo.url = url;
        geoInfo.error = ME.message;
    end
end

function structured = extractStructuredInfo(text)
    % Extract structured information using regex patterns
    structured = struct();
    
    % Define patterns for key fields
    patterns = struct();
    patterns.title = 'Title\s*\|\s*([^\|]+)';
    patterns.organism = 'Organism\s*\|\s*([^\|]+)';
    patterns.sample_type = 'Sample type\s*\|\s*([^\|]+)';
    patterns.source_name = 'Source name\s*\|\s*([^\|]+)';
    patterns.library_strategy = 'Library strategy\s*\|\s*([^\|]+)';
    patterns.instrument = 'Instrument model\s*\|\s*([^\|]+)';
    patterns.contact_name = 'Contact name\s*\|\s*([^\|]+)';
    patterns.organization = 'Organization name\s*\|\s*([^\|]+)';
    patterns.submission_date = 'Submission date\s*\|\s*([^\|]+)';
    
    % Extract information using patterns
    fields = fieldnames(patterns);
    for i = 1:length(fields)
        field = fields{i};
        pattern = patterns.(field);
        tokens = regexp(text, pattern, 'tokens', 'ignorecase');
        if ~isempty(tokens)
            structured.(field) = strtrim(tokens{1}{1});
        else
            structured.(field) = 'Not found';
        end
    end
    
    % Extract characteristics (more complex pattern)
    charPattern = 'Characteristics\s*\|\s*([^|]+?)(?=\s*\||$)';
    charTokens = regexp(text, charPattern, 'tokens', 'ignorecase');
    if ~isempty(charTokens)
        structured.characteristics = strtrim(charTokens{1}{1});
    else
        structured.characteristics = 'Not found';
    end
end

function llmAnalysis = analyzWithLLM(text, accession)
    % Use LLM to analyze the GEO record
    llmAnalysis = struct();
    
    try
        % Prepare prompt for LLM analysis
        prompt = sprintf([
            'Analyze this GEO database record (accession: %s) and extract key information:\n\n'...
            '%s\n\n'...
            'Please provide a structured analysis including:\n'...
            '1. Study type and experimental design\n'...
            '2. Biological significance\n'...
            '3. Technical details (sequencing platform, protocols)\n'...
            '4. Sample characteristics\n'...
            '5. Potential applications or research areas\n'...
            'Format your response as clear, concise bullet points.'
        ], accession, text);
        
        % Call LLM function
        % response = llm(prompt, 'ModelName', 'gpt-4', 'MaxTokens', 1000);
        % assignin('base',"prompt",prompt)

        response = llm.callTAMUAIChat([], prompt);
        llmAnalysis.detailed_analysis = response;
        
        % Extract research domain
        domainPrompt = sprintf([
            'Based on this GEO record, identify the primary research domain '...
            '(e.g., cancer biology, neuroscience, developmental biology, etc.):\n\n%s'
        ], text);
        
        % domainResponse = llm(domainPrompt, 'ModelName', 'gpt-4', 'MaxTokens', 100);
        domainResponse = llm.callTAMUAIChat([], domainPrompt);
        llmAnalysis.research_domain = domainResponse;
        
        % Extract key findings or objectives
        objectivePrompt = sprintf([
            'What are the main research objectives or expected findings from this study? '...
            'Provide a brief summary:\n\n%s'
        ], text);
        
        % objectiveResponse = llm(objectivePrompt, 'ModelName', 'gpt-4', 'MaxTokens', 200);
        objectiveResponse = llm.callTAMUAIChat([], objectivePrompt);
        llmAnalysis.research_objectives = objectiveResponse;
        
    catch ME
        warning(ME.identifier, 'LLM analysis failed: %s', ME.message);
        llmAnalysis.detailed_analysis = 'LLM analysis unavailable';
        llmAnalysis.research_domain = 'Unable to determine';
        llmAnalysis.research_objectives = 'Unable to extract';
    end
end

function textAnalysis = performTextAnalysis(text)
    % Use Text Analytics Toolbox for analysis
    textAnalysis = struct();
    
    try
        % Create text data
        textData = preprocessText(text);
        
        % Extract keywords using TF-IDF
        documents = tokenizedDocument(textData);
        bag = bagOfWords(documents);
        
        % Calculate TF-IDF scores
        tfidfScores = tfidf(bag);
        
        % Get top keywords
        [~, idx] = sort(tfidfScores, 'descend');
        topWords = bag.Vocabulary(idx(1:min(20, length(idx))));
        topScores = tfidfScores(idx(1:min(20, length(idx))));
        
        textAnalysis.top_keywords = topWords;
        textAnalysis.keyword_scores = topScores;
        
        % Extract named entities (if available)
        try
            entities = extractNamedEntities(textData);
            textAnalysis.named_entities = entities;
        catch
            textAnalysis.named_entities = 'Named entity extraction not available';
        end
        
        % Basic text statistics
        textAnalysis.word_count = length(split(text));
        textAnalysis.sentence_count = length(split(text, '.'));
        
    catch ME
        warning(ME.identifier, 'Text analysis failed: %s', ME.message);
        textAnalysis.top_keywords = {};
        textAnalysis.keyword_scores = [];
        textAnalysis.error = ME.message;
    end
end

function summary = generateSummary(geoInfo)
    % Generate a comprehensive summary
    summary = struct();
    
    % Basic information
    summary.accession = geoInfo.accession;
    summary.experiment_type = 'Unknown';
    summary.organism = 'Unknown';
    summary.key_findings = 'Analysis pending';
    
    % Extract key information from structured data
    if isfield(geoInfo, 'structured')
        s = geoInfo.structured;
        if ~strcmp(s.organism, 'Not found')
            summary.organism = s.organism;
        end
        if ~strcmp(s.library_strategy, 'Not found')
            summary.experiment_type = s.library_strategy;
        end
        if ~strcmp(s.title, 'Not found')
            summary.title = s.title;
        end
    end
    
    % Add LLM insights
    if isfield(geoInfo, 'llm_analysis')
        summary.research_domain = geoInfo.llm_analysis.research_domain;
        summary.research_objectives = geoInfo.llm_analysis.research_objectives;
    end
    
    % Add top keywords
    if isfield(geoInfo, 'text_analysis') && isfield(geoInfo.text_analysis, 'top_keywords')
        summary.top_keywords = geoInfo.text_analysis.top_keywords(1:min(10, length(geoInfo.text_analysis.top_keywords)));
    end
end

function plainText = extractHTMLText(htmlContent)
    % Simple HTML text extraction (you can use the previous function here)
    % Remove HTML tags
    plainText = regexprep(htmlContent, '<[^>]*>', '');
    
    % Convert HTML entities
    plainText = regexprep(plainText, '&amp;', '&');
    plainText = regexprep(plainText, '&lt;', '<');
    plainText = regexprep(plainText, '&gt;', '>');
    plainText = regexprep(plainText, '&quot;', '"');
    plainText = regexprep(plainText, '&nbsp;', ' ');
    
    % Clean up whitespace
    plainText = regexprep(plainText, '\s+', ' ');
    plainText = strtrim(plainText);
end