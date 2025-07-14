function plainText = ai_extractHTMLText(htmlInput)
    % EXTRACTHTMLTEXT Extract plain text from HTML content
    %
    % SYNTAX:
    %   plainText = extractHTMLText(htmlInput)
    %
    % INPUT:
    %   htmlInput - HTML content as a string, character array, or filename
    %
    % OUTPUT:
    %   plainText - Extracted plain text as a string
    %
    % EXAMPLES:
    %   % From HTML string
    %   html = '<p>Hello <b>world</b>!</p>';
    %   text = extractHTMLText(html);
    %
    %   % From HTML file
    %   text = extractHTMLText('webpage.html');
    
    % Convert input to string if needed
    if ischar(htmlInput)
        htmlInput = string(htmlInput);
    end
    
    % Check if input is a filename
    if isfile(htmlInput)
        % Read HTML from file
        htmlContent = fileread(htmlInput);
        htmlContent = string(htmlContent);
    else
        % Use input as HTML content directly
        htmlContent = htmlInput;
    end
    
    % Remove script and style tags and their content
    htmlContent = regexprep(htmlContent, '<script[^>]*>.*?</script>', '', 'ignorecase');
    htmlContent = regexprep(htmlContent, '<style[^>]*>.*?</style>', '', 'ignorecase');
    
    % Remove HTML comments
    htmlContent = regexprep(htmlContent, '<!--.*?-->', '');
    
    % Convert common HTML entities to characters
    htmlContent = regexprep(htmlContent, '&amp;', '&');
    htmlContent = regexprep(htmlContent, '&lt;', '<');
    htmlContent = regexprep(htmlContent, '&gt;', '>');
    htmlContent = regexprep(htmlContent, '&quot;', '"');
    htmlContent = regexprep(htmlContent, '&apos;', '''');
    htmlContent = regexprep(htmlContent, '&nbsp;', ' ');
    htmlContent = regexprep(htmlContent, '&#39;', '''');
    
    % Convert block-level elements to add line breaks
    blockElements = {'</p>', '</div>', '</h1>', '</h2>', '</h3>', '</h4>', '</h5>', '</h6>', ...
                     '</li>', '</tr>', '</td>', '</th>', '</blockquote>', '</pre>', '</br>', '<br>'};
    
    for i = 1:length(blockElements)
        htmlContent = regexprep(htmlContent, blockElements{i}, [blockElements{i} newline], 'ignorecase');
    end
    
    % Remove all remaining HTML tags
    htmlContent = regexprep(htmlContent, '<[^>]*>', '');
    
    % Clean up whitespace
    htmlContent = regexprep(htmlContent, '\s+', ' ');  % Multiple whitespace to single space
    htmlContent = regexprep(htmlContent, '^\s+|\s+$', '');  % Trim leading/trailing whitespace
    
    % Convert multiple newlines to single newlines
    htmlContent = regexprep(htmlContent, '\n\s*\n', '\n');
    
    % Remove empty lines
    lines = splitlines(htmlContent);
    lines = lines(~cellfun(@isempty, strtrim(lines)));
    
    % Join lines back together
    plainText = strjoin(lines, newline);
    
    % Final cleanup
    plainText = strtrim(plainText);
end