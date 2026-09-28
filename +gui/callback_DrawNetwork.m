function callback_DrawNetwork(src, ~)

[FigureHandle] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, "Input edge list from:","Select Source", ...
{'Paste Text', 'Open File', 'Cancel'},'Paste Text');
switch answer
    case 'Paste Text'
        defaulttxt = sprintf('ANGPT1\tITGB1\nRAC1\tVEGFA\nLAMC1\tDAG1\nLAMB1\tITGB1\nRAC1\tEZR\nMDK\tNOTCH2\nNID1\tPTPRF\nNID1\tITGB1\nLAMB1\tDAG1\nLAMC1\tITGB1\nTHBS1\tITGA3\nCOL6A3\tITGA3\nLAMA4\tITGA3\nCALR\tITGA3\nDCN\tVEGFA\nMDK\tVEGFA\nCALR\tPDIA3\nLAMA4\tDAG1\nRAC1\tITGB1\nCALR\tP4HB\nLAMB2\tITGB1\nSPON2\tITGB1\nVEGFA\tITGB1\nTHBS1\tCD47\nTLN1\tITGB1\nCOL16A1\tITGB1\nMDK\tNCL\n');
        % defaulttxt = sprintf('ANGPT1\tITGB1\nHSPA8\tADRB2\nRAC1\tVEGFA\nLAMC1\tDAG1\nLAMB1\tITGB1\nRAC1\tEZR\nMDK\tNOTCH2\nNID1\tPTPRF\nNID1\tITGB1\nLAMB1\tDAG1\nLAMC1\tITGB1\nTHBS1\tITGA3\nCOL6A3\tITGA3\nLAMA4\tITGA3\nCALR\tITGA3\nDCN\tVEGFA\nMDK\tVEGFA\nCALR\tPDIA3\nAPOE\tLRP6\nCXCL12\tGNAI2\nCTGF\tLRP6\nLAMA4\tDAG1\nCTGF\tF2RL1\nCD99\tCD81\nRAC1\tITGB1\nCALR\tP4HB\nLAMB2\tITGB1\nSPON2\tITGB1\nVEGFA\tITGB1\nTHBS1\tCD47\nSLIT2\tAPP\nTLN1\tITGB1\nVEGFB\tNRP1\nCOL16A1\tITGB1\nSEMA3A\tNRP1\nMDK\tNCL\n');
        % defaulttxt = sprintf('GeneA\tGeneB\nGeneA\tGeneC\nGeneA\tGeneD\nGeneB\tGeneD\n');

        if gui.i_isuifig(FigureHandle)
            % prompts = {sprintf(['Paste edge list\n'...
            %     'Format: Gene 1 [TAB] Gene 2'])};
            % [userInput] = gui.myInputdlg({sprintf(['Paste edge list\n' ...
            %     'Format: Gene 1 [TAB] Gene 2'])}, '', ...
            %     {defaulttxt}, FigureHandle);
            [userInput] = gui.myTextareadlg(FigureHandle, {''}, '', {defaulttxt}, true);
        else
            [userInput] = inputdlg(sprintf(['Paste edge list\n' ...
                'Format: Gene 1 [TAB] Gene 2']), '', ...
                [15, 80], {defaulttxt});
        end


        if isempty(userInput)
             disp('User canceled input.')
             return;
        end
        fw = gui.myWaitbar(FigureHandle);
        a = tempname;
        fid = fopen(a, 'w');
        fprintf(fid, '%s\n', string(userInput{1}));  % Write first input only (modify for multiple)
        fclose(fid);
        tab = readtable(a,"FileType","text", ...
            'Delimiter','\t','ReadVariableNames',false, ...
            'VariableNamingRule', 'modify');
    case 'Open File'
        % uigetfile takes no parent figure: a figure passed first is read as
        % the filter spec and errors. Raise the app again afterwards instead.
        [fname, pathname] = uigetfile( ...
            {'*.txt;*.tsv;*.mat', 'Network Files (*.txt, *.tsv, *.mat)'; ...
            '*.txt;*.tsv', 'Edge List (*.txt, *.tsv)'; ...
            '*.mat', 'Saved GRN (*.mat)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Pick a Network File');
        if pkg.i_isvalid(FigureHandle), figure(FigureHandle); end
        if isequal(fname, 0), return; end
        tabfile = fullfile(pathname, fname);
        [~, ~, ext] = fileparts(tabfile);
        if strcmpi(ext, '.mat')
            try
                [A, g, figname] = i_loadgrnmat(tabfile, FigureHandle);
            catch ME
                gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                return;
            end
            if isempty(A), return; end
            i_drawadjacency(A, g, figname, FigureHandle);
            return;
        end
        % Capture and restore rather than a bare off/on pair: the pair leaves
        % warnings disabled for the rest of the session if anything between the
        % two lines throws, and its 'on' re-enables warnings the caller may have
        % silenced deliberately instead of restoring what they had.
        warnState = warning();
        restoreWarn = onCleanup(@() warning(warnState));
        warning('off', 'all');
        fw = gui.myWaitbar(FigureHandle);
        tab = readtable(tabfile,'FileType','text', ...
            'Delimiter','\t','ReadVariableNames',false);
    otherwise
        return;
end
if width(tab) < 2
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ['The edge list needs two tab-separated ', ...
        'columns, Gene 1 [TAB] Gene 2, and optionally a third column of weights.']);
    return;
end
% A numeric third column is taken as the edge weight
if width(tab) >= 3 && isnumeric(tab.Var3)
    weights = double(tab.Var3);
    weights(isnan(weights)) = 1;
else
    weights = ones(height(tab), 1);
end
G = digraph(string(tab.Var1), string(tab.Var2), weights);
gui.myWaitbar(FigureHandle, fw);
if G.numnodes <= i_maxwholegenes()
    gui.i_singlegraph(G, '', FigureHandle);
else
    i_drawadjacency(adjacency(G, 'weighted'), string(G.Nodes.Name), '', ...
        FigureHandle);
end
end


function n = i_maxwholegenes()
% Largest network drawn gene by gene without asking; gui.i_showgrn uses the
% same limit for the networks it draws after a build
n = 100;
end


function i_drawadjacency(A, g, figname, FigureHandle)
% Draw a small network whole. A larger one is drawn as a map of its gene
% modules (sc_grnsketch), where a module or a single gene opens as a
% gene-level view; drawing it whole is offered only while that is still
% legible.
maxWholeEdges = 2000;
g = string(g(:));
numGenes = numel(g);
numEdges = nnz(A - diag(diag(A)));
if numGenes <= i_maxwholegenes()
    sc_grnview(A, g, figname, FigureHandle);
    return;
end
answer = 'Module Map';
if numEdges <= maxWholeEdges
    answer = gui.myQuestdlg(FigureHandle, sprintf(['The network has %d genes ', ...
        'and %d edges. Draw it as a map of gene modules (click a module or ', ...
        'find a gene to see its genes), or draw the whole network?'], ...
        numGenes, numEdges), 'Large Network', ...
        {'Module Map', 'Whole Network', 'Cancel'}, 'Module Map');
end
switch answer
    case 'Module Map'
        fw = gui.myWaitbar(FigureHandle, [], false, ...
            sprintf('Finding modules among %d genes...', numGenes));
        try
            sc_grnsketch(A, g, ParentFig=FigureHandle);
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        gui.myWaitbar(FigureHandle, fw);
    case 'Whole Network'
        sc_grnview(A, g, figname, FigureHandle);
    otherwise
        % Cancel or dialog closed
end
end


function [A, g, figname] = i_loadgrnmat(fname, FigureHandle)
% Read a network saved by scgeatool: A and g from Network > Build GRN,
% or A0, A1 and glist from Network > Compare GRNs. Also accepts genes or
% genelist for the gene list, and a graph or digraph variable. Returns an
% empty A when the user cancels. Only the needed variables are loaded: a
% comparison file also holds its result table.
A = [];
g = string.empty;
[~, figname] = fileparts(fname);
vars = string({whos('-file', fname).name});
graphVars = vars(arrayfun(@(v) any(strcmp(whos('-file', fname, v).class, ...
    {'graph', 'digraph'})), vars));

if any(vars == "A")
    varA = "A";
elseif all(ismember(["A0", "A1"], vars))
    groups = ["A0", "A1"];
    if any(vars == "groups")
        S = load(fname, 'groups');
        groups = string(S.groups(:)).';
    end
    answer = gui.myQuestdlg(FigureHandle, 'The file holds two networks. Draw which?', ...
        'Select Network', {char(groups(1)), char(groups(2)), 'Cancel'}, char(groups(1)));
    switch answer
        case char(groups(1))
            varA = "A0";
        case char(groups(2))
            varA = "A1";
        otherwise
            return;
    end
    figname = sprintf('%s (%s)', figname, answer);
elseif ~isempty(graphVars)
    S = load(fname, graphVars(1));
    G = S.(graphVars(1));
    A = adjacency(G, 'weighted');
    g = string(G.Nodes.Name);
    return;
else
    error('gui:callback_DrawNetwork:noNetwork', ['%s holds no network. ', ...
        'Expected an adjacency matrix A with gene list g, as saved by ', ...
        'Network > Build Gene Regulatory Network.'], fname);
end

geneVar = intersect(["g", "genes", "genelist", "glist"], vars, 'stable');
S = load(fname, varA, geneVar{:});
A = S.(varA);
if ~(isnumeric(A) || islogical(A)) || ~ismatrix(A) || size(A, 1) ~= size(A, 2)
    error('gui:callback_DrawNetwork:badAdjacency', ...
        '%s in %s is not a square adjacency matrix.', varA, fname);
end
if isempty(geneVar)
    g = "G" + string(1:size(A, 1)).';
else
    g = string(S.(geneVar(1)));
    g = g(:);
end
if numel(g) ~= size(A, 1)
    error('gui:callback_DrawNetwork:geneCount', ...
        'The gene list has %d entries but %s is %dx%d.', numel(g), varA, ...
        size(A, 1), size(A, 2));
end
end
