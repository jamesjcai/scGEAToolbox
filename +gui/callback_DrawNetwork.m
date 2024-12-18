function callback_DrawNetwork(src, ~)


import mlreportgen.ppt.*;

FigureHandle = src.Parent.Parent;

answer = questdlg("Input network table from:","Select Source", ...
    'Paste Text', 'Open File', 'Cancel','Paste Text');
switch answer
    case 'Paste Text'
    defaulttxt = sprintf('GeneA\tGeneB\nGeneA\tGeneC\nGeneA\tGeneD\nGeneB\tGeneD\n');
        [userInput] = inputdlg('Paste table text', 'Network', ...
            [15, 80], {defaulttxt});
        if isempty(userInput)
             disp('User canceled input.')
             return;
        end
        a = tempname;
        fid = fopen(a, 'w');
        fprintf(fid, '%s\n', string(userInput{1}));  % Write first input only (modify for multiple)
        fclose(fid);
        warning off
        tab = readtable(a,"FileType","text", ...
            'Delimiter','\t','ReadVariableNames',false);                
        warning on
    case 'Open File'
        [fname, pathname] = uigetfile( ...
            {'*.txt', 'Network Files (*.txt)'; ...
            '*.*', 'All Files (*.*)'}, ...
            'Pick a Network Table File');
        if isequal(fname, 0), return; end
        tabfile = fullfile(pathname, fname);
        warning off
        tab = readtable(tabfile,'FileType','text', ...
            'Delimiter','\t','ReadVariableNames',false);                
        warning on
    otherwise
        return;
end
weights = ones(height(tab), 1);
G = graph(tab.Var1, tab.Var2, weights);
gui.i_singlegraph(G, '', FigureHandle);
end