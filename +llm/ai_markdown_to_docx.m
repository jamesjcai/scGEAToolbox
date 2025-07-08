function ai_markdown_to_docx
    % MATLAB GUI for converting Markdown to DOCX/PDF with improved functionality
    import mlreportgen.dom.*

    % Create main figure
    fig = uifigure('Name', 'Markdown to DOCX & PDF Converter', ...
                   'Position', [100, 100, 700, 650], ...
                   'Resize', 'on');

    % Title label
    titleLabel = uilabel(fig, 'Position', [30, 600, 640, 30], ...
                        'Text', 'Markdown to DOCX & PDF Converter', ...
                        'FontSize', 16, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center');

    % Input section
    inputLabel = uilabel(fig, 'Position', [30, 570, 250, 20], ...
                        'Text', 'Enter your Markdown content:');
    
    % Markdown input area with improved default content
    markdownInput = uitextarea(fig, 'Position', [30, 320, 640, 240], ...
                              'FontName', 'Consolas');
    markdownInput.Value = {
        '# Sample Document',
        '',
        'This is a **bold** statement and this is *italic* text.',
        'You can also have ***bold and italic*** text.',
        '',
        '## Features',
        '',
        '- Supports headings (H1-H6)',
        '- Text formatting: **bold**, *italic*, ***bold+italic***',
        '- Ordered and unordered lists',
        '- Multiple paragraphs',
        '',
        '### Ordered List Example',
        '',
        '1. First item',
        '2. Second item',
        '3. Third item',
        '',
        'Regular paragraph text continues here.'
    };

    % Preview section
    previewLabel = uilabel(fig, 'Position', [30, 290, 250, 20], ...
                          'Text', 'Preview:');
    
    % Preview area
    previewArea = uihtml(fig, 'Position', [30, 90, 640, 190], ...
                        'HTMLSource', '');

    % Buttons
    btnPreview = uibutton(fig, 'Text', 'Preview', ...
                         'Position', [30, 40, 100, 35], ...
                         'ButtonPushedFcn', @(~,~) preview_markdown());
    
    btnClear = uibutton(fig, 'Text', 'Clear', ...
                       'Position', [150, 40, 100, 35], ...
                       'ButtonPushedFcn', @(~,~) clear_content());
    
    btnDOCX = uibutton(fig, 'Text', 'Export DOCX', ...
                      'Position', [380, 40, 120, 35], ...
                      'ButtonPushedFcn', @(~,~) export_file('docx'));
    
    btnPDF = uibutton(fig, 'Text', 'Export PDF', ...
                     'Position', [520, 40, 120, 35], ...
                     'ButtonPushedFcn', @(~,~) export_file('pdf'));

    % Status label
    statusLabel = uilabel(fig, 'Position', [30, 10, 640, 20], ...
                         'Text', 'Ready', 'FontColor', [0, 0.5, 0]);

    % Initialize with preview
    preview_markdown();

    function preview_markdown()
        try
            statusLabel.Text = 'Generating preview...';
            statusLabel.FontColor = [0, 0, 1];
            drawnow;
            
            md = markdownInput.Value;
            if ischar(md)
                md = {md};
            end
            
            htmlOutput = convert_markdown_to_html(md);
            
            % Create styled HTML with better formatting
            styledHTML = sprintf(['<html><head><style>', ...
                'body { font-family: Arial, sans-serif; font-size: 14px; ', ...
                'line-height: 1.6; color: black; padding: 10px; }', ...
                'h1 { font-size: 24px; margin: 20px 0 10px 0; }', ...
                'h2 { font-size: 20px; margin: 18px 0 8px 0; }', ...
                'h3 { font-size: 16px; margin: 16px 0 6px 0; }', ...
                'h4, h5, h6 { font-size: 14px; margin: 14px 0 4px 0; }', ...
                'p { margin: 10px 0; }', ...
                'ul, ol { margin: 10px 0; padding-left: 30px; }', ...
                'li { margin: 4px 0; }', ...
                '</style></head><body>%s</body></html>'], htmlOutput);
            
            previewArea.HTMLSource = styledHTML;
            
            statusLabel.Text = 'Preview updated successfully';
            statusLabel.FontColor = [0, 0.5, 0];
            
        catch ME
            statusLabel.Text = sprintf('Preview error: %s', ME.message);
            statusLabel.FontColor = [1, 0, 0];
        end
    end

    function clear_content()
        markdownInput.Value = {''};
        previewArea.HTMLSource = '';
        statusLabel.Text = 'Content cleared';
        statusLabel.FontColor = [0, 0.5, 0];
    end

    function export_file(format)
        try
            statusLabel.Text = sprintf('Exporting %s...', upper(format));
            statusLabel.FontColor = [0, 0, 1];
            drawnow;
            
            md = markdownInput.Value;
            if ischar(md)
                md = {md};
            end
            
            % Get save location
            [filename, pathname] = uiputfile(['*.', format], ...
                                           ['Save as ', upper(format)]);
            if isequal(filename, 0)
                statusLabel.Text = 'Export cancelled';
                statusLabel.FontColor = [0.5, 0.5, 0.5];
                return;
            end
            
            fullPath = fullfile(pathname, filename);
            
            % Create document
            doc = Document(fullPath, upper(format));
            
            % Parse markdown and add to document
            domBlocks = parse_markdown_to_dom(md);
            
            for i = 1:length(domBlocks)
                append(doc, domBlocks{i});
            end
            
            close(doc);
            
            % Open the file
            if strcmpi(format, 'pdf')
                web(fullPath, '-browser');
            else
                winopen(fullPath);
            end
            
            statusLabel.Text = sprintf('%s exported successfully: %s', ...
                                     upper(format), filename);
            statusLabel.FontColor = [0, 0.5, 0];
            
        catch ME
            statusLabel.Text = sprintf('Export error: %s', ME.message);
            statusLabel.FontColor = [1, 0, 0];
        end
    end
end

function domBlocks = parse_markdown_to_dom(markdownLines)
    % Parse markdown content into DOM blocks
    import mlreportgen.dom.*
    
    if ischar(markdownLines)
        markdownLines = {markdownLines};
    end
    
    % Join lines and split again for consistent processing
    content = strjoin(markdownLines, newline);
    lines = strsplit(content, newline);
    
    domBlocks = {};
    i = 1;
    
    while i <= length(lines)
        line = strtrim(lines{i});
        
        % Skip empty lines
        if isempty(line)
            i = i + 1;
            continue;
        end
        
        % Check for headers
        headerMatch = regexp(line, '^(#{1,6})\s+(.*)$', 'tokens', 'once');
        if ~isempty(headerMatch)
            level = length(headerMatch{1});
            level = min(level, 6); % Ensure level doesn't exceed 6
            content = parse_inline_formatting(headerMatch{2});
            domBlocks{end+1} = Heading(level, content);
            i = i + 1;
            continue;
        end
        
        % Check for ordered list
        olMatch = regexp(line, '^\d+\.\s+(.*)$', 'tokens', 'once');
        if ~isempty(olMatch)
            ol = OrderedList();
            while i <= length(lines)
                currentLine = strtrim(lines{i});
                currentMatch = regexp(currentLine, '^\d+\.\s+(.*)$', 'tokens', 'once');
                if isempty(currentMatch)
                    break;
                end
                content = parse_inline_formatting(currentMatch{1});
                append(ol, ListItem(content));
                i = i + 1;
            end
            domBlocks{end+1} = ol;
            continue;
        end
        
        % Check for unordered list
        ulMatch = regexp(line, '^[-*+]\s+(.*)$', 'tokens', 'once');
        if ~isempty(ulMatch)
            ul = UnorderedList();
            while i <= length(lines)
                currentLine = strtrim(lines{i});
                currentMatch = regexp(currentLine, '^[-*+]\s+(.*)$', 'tokens', 'once');
                if isempty(currentMatch)
                    break;
                end
                content = parse_inline_formatting(currentMatch{1});
                append(ul, ListItem(content));
                i = i + 1;
            end
            domBlocks{end+1} = ul;
            continue;
        end
        
        % Regular paragraph
        content = parse_inline_formatting(line);
        domBlocks{end+1} = Paragraph(content);
        i = i + 1;
    end
end

function content = parse_inline_formatting(text)
    % Parse inline formatting (bold, italic, bold+italic)
    import mlreportgen.dom.*
    
    if isempty(text)
        content = Text('');
        return;
    end
    
    % Pattern to match formatting: ***text***, **text**, *text*
    pattern = '(\*\*\*[^*]+\*\*\*|\*\*[^*]+\*\*|\*[^*]+\*)';
    
    % Find all matches
    [startIdx, endIdx, matches] = regexp(text, pattern, 'start', 'end', 'match');
    
    if isempty(matches)
        % No formatting found
        content = Text(text);
        return;
    end
    
    % Build content array with formatted and unformatted parts
    content = {};
    lastEnd = 0;
    
    for i = 1:length(matches)
        % Add text before the match
        if startIdx(i) > lastEnd + 1
            content{end+1} = Text(text(lastEnd+1:startIdx(i)-1));
        end
        
        % Process the match
        match = matches{i};
        if startsWith(match, '***') && endsWith(match, '***')
            % Bold and italic
            innerText = match(4:end-3);
            t = Text(innerText);
            t.Bold = true;
            t.Italic = true;
            content{end+1} = t;
        elseif startsWith(match, '**') && endsWith(match, '**')
            % Bold only
            innerText = match(3:end-2);
            t = Text(innerText);
            t.Bold = true;
            content{end+1} = t;
        elseif startsWith(match, '*') && endsWith(match, '*')
            % Italic only
            innerText = match(2:end-1);
            t = Text(innerText);
            t.Italic = true;
            content{end+1} = t;
        end
        
        lastEnd = endIdx(i);
    end
    
    % Add remaining text
    if lastEnd < length(text)
        content{end+1} = Text(text(lastEnd+1:end));
    end
    
    % Return single element if only one, otherwise return array
    if length(content) == 1
        content = content{1};
    end
end

function html = convert_markdown_to_html(markdownLines)
    % Convert markdown to HTML for preview
    if ischar(markdownLines)
        markdownLines = {markdownLines};
    end
    
    html = '';
    i = 1;
    
    while i <= length(markdownLines)
        line = strtrim(markdownLines{i});
        
        if isempty(line)
            html = [html, '<br>'];
            i = i + 1;
            continue;
        end
        
        % Headers
        headerMatch = regexp(line, '^(#{1,6})\s+(.*)$', 'tokens', 'once');
        if ~isempty(headerMatch)
            level = length(headerMatch{1});
            content = format_inline_html(headerMatch{2});
            html = [html, sprintf('<h%d>%s</h%d>', level, content, level)];
            i = i + 1;
            continue;
        end
        
        % Ordered list
        olMatch = regexp(line, '^\d+\.\s+(.*)$', 'tokens', 'once');
        if ~isempty(olMatch)
            html = [html, '<ol>'];
            while i <= length(markdownLines)
                currentLine = strtrim(markdownLines{i});
                currentMatch = regexp(currentLine, '^\d+\.\s+(.*)$', 'tokens', 'once');
                if isempty(currentMatch)
                    break;
                end
                content = format_inline_html(currentMatch{1});
                html = [html, sprintf('<li>%s</li>', content)];
                i = i + 1;
            end
            html = [html, '</ol>'];
            continue;
        end
        
        % Unordered list
        ulMatch = regexp(line, '^[-*+]\s+(.*)$', 'tokens', 'once');
        if ~isempty(ulMatch)
            html = [html, '<ul>'];
            while i <= length(markdownLines)
                currentLine = strtrim(markdownLines{i});
                currentMatch = regexp(currentLine, '^[-*+]\s+(.*)$', 'tokens', 'once');
                if isempty(currentMatch)
                    break;
                end
                content = format_inline_html(currentMatch{1});
                html = [html, sprintf('<li>%s</li>', content)];
                i = i + 1;
            end
            html = [html, '</ul>'];
            continue;
        end
        
        % Regular paragraph
        content = format_inline_html(line);
        html = [html, sprintf('<p>%s</p>', content)];
        i = i + 1;
    end
end

function formatted = format_inline_html(text)
    % Format inline markdown to HTML
    if isempty(text)
        formatted = '';
        return;
    end
    
    formatted = text;
    
    % Bold and italic: ***text***
    formatted = regexprep(formatted, '\*\*\*([^*]+)\*\*\*', '<strong><em>$1</em></strong>');
    
    % Bold: **text**
    formatted = regexprep(formatted, '\*\*([^*]+)\*\*', '<strong>$1</strong>');
    
    % Italic: *text*
    formatted = regexprep(formatted, '\*([^*]+)\*', '<em>$1</em>');
end