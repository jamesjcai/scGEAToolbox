function plainText = extractplaintext(htmlCode)
    % EXTRACTPLAINTEXT Extract plain text from HTML code
    %
    % Syntax:
    %   plainText = extractplaintext(htmlCode)
    %
    % Input:
    %   htmlCode - string or char array containing HTML code
    %
    % Output:
    %   plainText - string containing extracted plain text
    %
    % Example:
    %   html = '<html><body><h1>Title</h1><p>This is a paragraph.</p></body></html>';
    %   text = extractplaintext(html);
    
    % Convert to string if input is char array
    if ischar(htmlCode)
        htmlCode = string(htmlCode);
    end
    
    % Remove script and style tags and their content
    htmlCode = regexprep(htmlCode, '<script[^>]*>.*?</script>', '', 'ignorecase');
    htmlCode = regexprep(htmlCode, '<style[^>]*>.*?</style>', '', 'ignorecase');
    
    % Remove HTML comments
    htmlCode = regexprep(htmlCode, '<!--.*?-->', '');
    
    % Replace common block elements with newlines to preserve structure
    blockElements = {'</div>', '</p>', '</h1>', '</h2>', '</h3>', '</h4>', '</h5>', '</h6>', ...
                     '</li>', '</tr>', '</td>', '</th>', '<br>', '<br/>', '<br />'};
    
    for i = 1:length(blockElements)
        htmlCode = regexprep(htmlCode, blockElements{i}, '\n', 'ignorecase');
    end
    
    % Remove all remaining HTML tags
    htmlCode = regexprep(htmlCode, '<[^>]*>', '');
    
    % Decode common HTML entities
    htmlEntities = containers.Map({...
        '&amp;', '&lt;', '&gt;', '&quot;', '&apos;', '&nbsp;', ...
        '&copy;', '&reg;', '&trade;', '&hellip;', '&mdash;', '&ndash;', ...
        '&ldquo;', '&rdquo;', '&lsquo;', '&rsquo;'}, {...
        '&', '<', '>', '"', '''', ' ', ...
        '©', '®', '™', '...', '—', '–', ...
        '"', '"', '''', ''''});
    
    entityKeys = keys(htmlEntities);
    for i = 1:length(entityKeys)
        htmlCode = regexprep(htmlCode, entityKeys{i}, htmlEntities(entityKeys{i}), 'ignorecase');
    end
    
    % Decode numeric HTML entities (e.g., &#39; &#x27;)
    htmlCode = regexprep(htmlCode, '&#(\d+);', '${char(str2double($1))}');
    htmlCode = regexprep(htmlCode, '&#x([0-9a-fA-F]+);', '${char(hex2dec($1))}');
    
    % Clean up whitespace
    % Replace multiple consecutive whitespace characters with single space
    htmlCode = regexprep(htmlCode, '\s+', ' ');
    
    % Replace multiple consecutive newlines with double newline
    htmlCode = regexprep(htmlCode, '\n\s*\n\s*\n+', '\n\n');
    
    % Trim leading and trailing whitespace
    plainText = strip(htmlCode);
    
    % Convert back to string type
    plainText = string(plainText);
end