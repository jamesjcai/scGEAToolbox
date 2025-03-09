function callback_DrawNetwork(src, ~)


import mlreportgen.ppt.*;

[FigureHandle, ~, isui] = gui.gui_getfigsce(src);

answer = questdlg("Input edge list from:","Select Source", ...
    'Paste Text', 'Open File', 'Cancel','Paste Text');
switch answer
    case 'Paste Text'
        defaulttxt = sprintf('ANGPT1\tITGB1\nRAC1\tVEGFA\nLAMC1\tDAG1\nLAMB1\tITGB1\nRAC1\tEZR\nMDK\tNOTCH2\nNID1\tPTPRF\nNID1\tITGB1\nLAMB1\tDAG1\nLAMC1\tITGB1\nTHBS1\tITGA3\nCOL6A3\tITGA3\nLAMA4\tITGA3\nCALR\tITGA3\nDCN\tVEGFA\nMDK\tVEGFA\nCALR\tPDIA3\nLAMA4\tDAG1\nRAC1\tITGB1\nCALR\tP4HB\nLAMB2\tITGB1\nSPON2\tITGB1\nVEGFA\tITGB1\nTHBS1\tCD47\nTLN1\tITGB1\nCOL16A1\tITGB1\nMDK\tNCL\n');
        %defaulttxt = sprintf('ANGPT1\tITGB1\nHSPA8\tADRB2\nRAC1\tVEGFA\nLAMC1\tDAG1\nLAMB1\tITGB1\nRAC1\tEZR\nMDK\tNOTCH2\nNID1\tPTPRF\nNID1\tITGB1\nLAMB1\tDAG1\nLAMC1\tITGB1\nTHBS1\tITGA3\nCOL6A3\tITGA3\nLAMA4\tITGA3\nCALR\tITGA3\nDCN\tVEGFA\nMDK\tVEGFA\nCALR\tPDIA3\nAPOE\tLRP6\nCXCL12\tGNAI2\nCTGF\tLRP6\nLAMA4\tDAG1\nCTGF\tF2RL1\nCD99\tCD81\nRAC1\tITGB1\nCALR\tP4HB\nLAMB2\tITGB1\nSPON2\tITGB1\nVEGFA\tITGB1\nTHBS1\tCD47\nSLIT2\tAPP\nTLN1\tITGB1\nVEGFB\tNRP1\nCOL16A1\tITGB1\nSEMA3A\tNRP1\nMDK\tNCL\n');
        %defaulttxt = sprintf('GeneA\tGeneB\nGeneA\tGeneC\nGeneA\tGeneD\nGeneB\tGeneD\n');
        [userInput] = inputdlg(sprintf('Paste edge list\nFormat: Gene 1 [TAB] Gene 2'), '', ...
            [15, 80], {defaulttxt});
        if isempty(userInput)
             disp('User canceled input.')
             return;
        end
        fw = gui.gui_waitbar;
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
        fw = gui.gui_waitbar;
        tab = readtable(tabfile,'FileType','text', ...
            'Delimiter','\t','ReadVariableNames',false);                
        warning on
    otherwise
        return;
end
weights = ones(height(tab), 1);
G = digraph(tab.Var1, tab.Var2, weights);
gui.i_singlegraph(G, '', FigureHandle);
gui.gui_waitbar(fw);
end