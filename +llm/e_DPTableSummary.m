function [done, outfile] = e_DPTableSummary(Tup, Tdn, infotagstr, wrkdir)
% E_DPTABLESUMMARY - Generate LLM summary report directly from DP analysis results
%
% Unlike DE/DV analyses, DP results are already at the pathway/program level,
% so no Enrichr step is needed. The LLM directly interprets the differential
% program activity table produced by sc_dpg.
%
% Inputs:
%   Tup         : table of up-regulated programs (avg_log2FC > 0) from sc_dpg
%   Tdn         : table of down-regulated programs (avg_log2FC < 0) from sc_dpg
%   infotagstr  : string describing the comparison (used in prompt + output filename)
%   wrkdir      : output directory for the Word document

done = false;
outfile = [];

if nargin < 3, infotagstr = pkg.i_randinfostr; end
if nargin < 4, wrkdir = []; end

if ~isfolder(wrkdir)
    wrkdir = getpref('scgeatoolbox', 'externalwrkpath');
    if isempty(wrkdir), return; end
end

preftagname = 'llmodelprovider';
s = getpref('scgeatoolbox', preftagname);
providermodel = strsplit(s, ':');

% Build JSON for up- and down-regulated programs (program name, fold-change, p-value)
s_up = '';
s_dn = '';
if ~isempty(Tup)
    Tout = Tup(:, {'setnames', 'avg_log2FC', 'p_val_adj'});
    Tout.Properties.VariableNames = {'Program', 'log2FC', 'adj_pval'};
    s_up = jsonencode(table2struct(Tout));
end
if ~isempty(Tdn)
    Tout = Tdn(:, {'setnames', 'avg_log2FC', 'p_val_adj'});
    Tout.Properties.VariableNames = {'Program', 'log2FC', 'adj_pval'};
    s_dn = jsonencode(table2struct(Tout));
end

if isempty(s_up) && isempty(s_dn), return; end

prompt1 = "You are a researcher analyzing differential gene program activity in single-cell RNA-seq data. " + ...
    "Gene programs are predefined sets of biologically related genes, such as signaling pathway response modules or transcription factor targets. " + ...
    "You will be given a list of gene programs that are significantly more or less active in one cell group compared to another, " + ...
    "along with their log2 fold-changes and adjusted p-values. " + ...
    "The comparison being analyzed is: " + infotagstr + ". " + ...
    "Please summarize the biological significance of these differential programs in plain text, " + ...
    "highlighting the key pathways and what they suggest about the biological differences between the two groups. " + ...
    "Please write an executive summary to report the results of your analysis. ";

feedbk_up = [];
feedbk_dn = [];
done1 = false;
done2 = false;

switch providermodel{1}
    case 'Ollama'
        try
            webread("http://localhost:11434");
        catch
            disp('Ollama is not running.');
            return;
        end
        if ~isempty(s_up)
            [done1, feedbk_up] = llm.callOllama(prompt1 + "Up-regulated programs: " + s_up, providermodel{2});
        end
        if ~isempty(s_dn)
            [done2, feedbk_dn] = llm.callOllama(prompt1 + "Down-regulated programs: " + s_dn, providermodel{2});
        end
    case 'TAMUAIChat'
        if ~isempty(s_up)
            [done1, feedbk_up] = llm.callTAMUAIChat([], prompt1 + "Up-regulated programs: " + s_up, providermodel{2});
        end
        if ~isempty(s_dn)
            [done2, feedbk_dn] = llm.callTAMUAIChat([], prompt1 + "Down-regulated programs: " + s_dn, providermodel{2});
        end
    case 'Gemini'
        if ~isempty(s_up)
            [done1, feedbk_up] = llm.callGemini([], prompt1 + "Up-regulated programs: " + s_up, providermodel{2});
        end
        if ~isempty(s_dn)
            [done2, feedbk_dn] = llm.callGemini([], prompt1 + "Down-regulated programs: " + s_dn, providermodel{2});
        end
    case 'OpenAI'
        if ~isempty(s_up)
            [done1, feedbk_up] = llm.callOpenAIChat([], prompt1 + "Up-regulated programs: " + s_up, providermodel{2});
        end
        if ~isempty(s_dn)
            [done2, feedbk_dn] = llm.callOpenAIChat([], prompt1 + "Down-regulated programs: " + s_dn, providermodel{2});
        end
    otherwise
        warning('Invalid LLM provider and/or model.');
        return;
end

if ~(done1 || done2)
    fprintf('LLM call failed for "%s". Check provider settings and API key.\n', infotagstr);
    return;
end
if isempty(feedbk_up) && isempty(feedbk_dn), return; end

import mlreportgen.dom.*
outfile = fullfile(wrkdir, "Res_" + matlab.lang.makeValidName(infotagstr));
doc = mlreportgen.dom.Document(outfile, 'docx');
open(doc);
i_todoc(doc, sprintf('%s — Up-regulated programs', infotagstr), feedbk_up);
i_todoc(doc, sprintf('%s — Down-regulated programs', infotagstr), feedbk_dn);
close(doc);
done = true;

function i_todoc(doc, titstr, text)
    if isempty(text), return; end
    import mlreportgen.dom.*
    p = Paragraph(titstr);
    p.Style = {Bold(true), FontSize('14pt'), Color('blue')};
    append(doc, p);
    try
        text = regexprep(text, '<think>[\s\S]*?</think>', '');
    catch ME
        disp(ME.message);
    end
    try
        domBlocks = llm.i_parsemarkdowntodom(text);
        for k = 1:numel(domBlocks)
            append(doc, domBlocks{k});
        end
    catch ME
        disp(ME.message);
    end
end

end
