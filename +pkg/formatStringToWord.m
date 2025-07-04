function formatStringToWord(inputString, filename, varargin)
% formatStringToWord - Formats a string and saves it to a Word document
%
% Syntax:
%   formatStringToWord(inputString, filename)
%   formatStringToWord(inputString, filename, 'PropertyName', PropertyValue, ...)
%
% Inputs:
%   inputString - The text string to format (supports markdown-like syntax)
%   filename    - Output filename (with or without .docx extension)
%
% Optional Parameters:
%   'FontFamily'     - Default font family (default: 'Arial')
%   'FontSize'       - Default font size (default: '11pt')
%   'HighlightGenes' - Whether to highlight gene names (default: true)
%   'AutoOpen'       - Whether to open document after creation (default: false)
%   'Title'          - Document title (default: extracted from first line)
%
% Example:
%   text = "# My Report\n## Section 1\nThis is body text with STAT1 gene.";
%   formatStringToWord(text, 'MyReport.docx', 'HighlightGenes', true);

% Import required classes
import mlreportgen.dom.*

% Parse input arguments
p = inputParser;
addRequired(p, 'inputString', @ischar);
addRequired(p, 'filename', @ischar);
addParameter(p, 'FontFamily', 'Arial', @ischar);
addParameter(p, 'FontSize', '11pt', @ischar);
addParameter(p, 'HighlightGenes', true, @islogical);
addParameter(p, 'AutoOpen', false, @islogical);
addParameter(p, 'Title', '', @ischar);

parse(p, inputString, filename, varargin{:});

% Extract parameters
fontFamily = p.Results.FontFamily;
fontSize = p.Results.FontSize;
highlightGenes = p.Results.HighlightGenes;
autoOpen = p.Results.AutoOpen;
docTitle = p.Results.Title;

% Ensure filename has .docx extension
[~, name, ext] = fileparts(filename);
if isempty(ext)
    filename = [filename '.docx'];
elseif ~strcmpi(ext, '.docx')
    filename = [name '.docx'];
end

% Create document
doc = Document(filename, 'docx');

% Define styles
titleStyle = {FontFamily(fontFamily), FontSize('18pt'), Bold, Color('blue')};
heading1Style = {FontFamily(fontFamily), FontSize('16pt'), Bold, Color('darkblue')};
heading2Style = {FontFamily(fontFamily), FontSize('14pt'), Bold, Color('navy')};
heading3Style = {FontFamily(fontFamily), FontSize('12pt'), Bold, Color('black')};
bodyStyle = {FontFamily(fontFamily), FontSize(fontSize), Color('black')};
geneStyle = {FontFamily(fontFamily), FontSize(fontSize), Bold, Color('darkgreen')};

% Split text into lines
lines = strsplit(inputString, '\n');

% Process each line
for i = 1:length(lines)
    line = strtrim(lines{i});
    
    if isempty(line)
        % Add empty paragraph for spacing
        append(doc, Paragraph(' '));
        continue;
    end
    
    % Determine line type and format accordingly
    if startsWith(line, '# ')
        % Main title (H1)
        title = line(3:end);
        if isempty(docTitle)
            docTitle = title;
        end
        para = Paragraph(title);
        para.Style = titleStyle;
        append(doc, para);
        
    elseif startsWith(line, '## ')
        % Heading 1 (H2)
        heading = line(4:end);
        para = Paragraph(heading);
        para.Style = heading1Style;
        append(doc, para);
        
    elseif startsWith(line, '### ')
        % Heading 2 (H3)
        heading = line(5:end);
        para = Paragraph(heading);
        para.Style = heading2Style;
        append(doc, para);
        
    elseif startsWith(line, '#### ')
        % Heading 3 (H4)
        heading = line(6:end);
        para = Paragraph(heading);
        para.Style = heading3Style;
        append(doc, para);
        
    else
        % Body text - process for gene highlighting
        if highlightGenes
            para = processBodyTextWithGenes(line, bodyStyle, geneStyle);
        else
            para = Paragraph(line);
            para.Style = bodyStyle;
        end
        append(doc, para);
    end
end

% Add footer with metadata
append(doc, Paragraph(' '));
append(doc, Paragraph(' '));
footerText = Paragraph(['Document generated on: ' datestr(now, 'yyyy-mm-dd HH:MM:SS')]);
footerText.Style = {FontFamily(fontFamily), FontSize('9pt'), Color('gray'), Italic};
append(doc, footerText);

% Close and save document
close(doc);

% Display completion message
fprintf('Word document "%s" has been created successfully.\n', filename);
fprintf('Location: %s\n', doc.OutputPath);

% Auto-open if requested
if autoOpen
    try
        winopen(doc.OutputPath);
    catch
        fprintf('Could not auto-open document. Please open manually.\n');
    end
end

end

function para = processBodyTextWithGenes(text, bodyStyle, geneStyle)
% Process body text to highlight gene names
import mlreportgen.dom.*

% Common gene patterns (add more as needed)
genePatterns = {
    'STAT[0-9]+', 'IFITM[0-9]+', 'ISG[0-9]+', 'SP[0-9]+', 'BST[0-9]+', ...
    'IFIT[0-9]+', 'ZBP[0-9]+', 'SAMHD[0-9]+', 'IRF[0-9]+', 'MX[0-9]+', ...
    'OAS[0-9]+', 'RSAD[0-9]+', 'TRIM[0-9]+', 'PARP[0-9]+', 'DDX[0-9]+', ...
    'HERC[0-9]+', 'USP[0-9]+', 'PLSCR[0-9]+', 'GBP[0-9]+', 'CXCL[0-9]+', ...
    'CCL[0-9]+', 'IL[0-9]+', 'TNF[A-Z]*', 'IFNG', 'IFNA', 'IFNB'
};

% Create combined pattern
combinedPattern = strjoin(genePatterns, '|');

% Find gene matches - use 'start' and 'end' for more reliable indexing
[matchStart, matchEnd] = regexp(text, combinedPattern, 'start', 'end');
matches = regexp(text, combinedPattern, 'match');

% Create paragraph
para = Paragraph();
para.Style = bodyStyle;

if isempty(matches)
    % No genes found, add as regular text
    append(para, Text(text));
else
    % Process text with gene highlighting
    lastEnd = 1;
    
    for i = 1:length(matches)
        startIdx = matchStart(i);
        endIdx = matchEnd(i);
        
        % Add text before gene
        if startIdx > lastEnd
            beforeText = text(lastEnd:startIdx-1);
            if ~isempty(beforeText)
                append(para, Text(beforeText));
            end
        end
        
        % Add highlighted gene
        geneText = Text(matches{i});
        geneText.Style = geneStyle;
        append(para, geneText);
        
        lastEnd = endIdx + 1;
    end
    
    % Add remaining text
    if lastEnd <= length(text)
        remainingText = text(lastEnd:end);
        if ~isempty(remainingText)
            append(para, Text(remainingText));
        end
    end
end

end

% Example usage function
function exampleUsage()
% Example of how to use the formatStringToWord function

% Example 1: Simple text with markdown-like headers
text1 = ['# Executive Summary: Gene Function Enrichment Analysis', newline, ...
         '## Overview', newline, ...
         'This analysis reveals immune response signatures with STAT1 and STAT2 genes.', newline, ...
         '### Key Findings', newline, ...
         'Important genes include IFITM1, ISG15, and SP100 in interferon pathways.'];

formatStringToWord(text1, 'example1.docx', 'HighlightGenes', true);

% Example 2: Your original text
text2 = ['# Executive Summary: Gene Function Enrichment Analysis', newline, ...
         '## Overview', newline, ...
         'This analysis of gene function enrichment results reveals a strong immune response signature, particularly centered around interferon signaling, antiviral defense mechanisms, and inflammatory processes. The gene set shows significant enrichment in biological pathways related to cytokine-mediated signaling, antimicrobial responses, and regulation of immune cell migration and proliferation.', newline, ...
         '## Key Biological Processes', newline, ...
         '### Interferon Response Pathways', newline, ...
         'The gene set shows substantial enrichment in interferon response pathways, particularly Type I and Type II interferon signaling. Key genes include STAT1, STAT2, IFITM1, ISG15, and SP100. These genes are crucial for cellular antiviral defense mechanisms and immune system regulation.', newline, ...
         '### Antiviral Defense Mechanisms', newline, ...
         'A prominent feature of the gene set is its strong association with viral defense responses. Genes such as BST2, IFIT1, IFIT3, ZBP1, and SAMHD1 are involved in defense response to viruses and negative regulation of viral processes.'];

formatStringToWord(text2, 'gene_analysis.docx', 'HighlightGenes', true, 'AutoOpen', false);

fprintf('Example documents created successfully!\n');
end

% Utility function to read text from file
function formatTextFileToWord(textFilePath, outputFilename, varargin)
% Read text from file and format to Word document
%
% Syntax:
%   formatTextFileToWord(textFilePath, outputFilename)
%   formatTextFileToWord(textFilePath, outputFilename, 'PropertyName', PropertyValue, ...)

% Read text from file
if exist(textFilePath, 'file')
    fid = fopen(textFilePath, 'r', 'n', 'UTF-8');
    inputString = fread(fid, '*char')';
    fclose(fid);
    
    % Format to Word
    formatStringToWord(inputString, outputFilename, varargin{:});
else
    error('Text file not found: %s', textFilePath);
end

end