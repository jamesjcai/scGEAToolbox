function callback_EnrichrTab2Circos(src, ~, tab)

if nargin<3, tab = []; end

import mlreportgen.ppt.*;
%pw1 = fileparts(mfilename('fullpath'));
%pth = fullfile(pw1, '..', 'assets', 'myTemplate.pptx');

[FigureHandle] = gui.gui_getfigsce(src);

if isempty(tab)
answer1 = gui.myQuestdlg(FigureHandle, "Select the source of Enrichr result table.","", ...
    {'Web Enrichr', 'Matlab Enrichr'}, 'Web Enrichr');

switch answer1 
    case 'Matlab Enrichr'
        a = evalin('base', 'whos');
        b = struct2cell(a);
        valididx = ismember(b(4, :), 'table');
        if sum(valididx) < 1
            gui.myWarndlg(FigureHandle, 'No table variables in Workspace.');
            return;
        end
        b = b(:, valididx);
        a = a(valididx);
        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, b(1, :), 'Select Enrichr Result Table:');
        else
            [indx, tf] = listdlg('PromptString', {'Select Enrichr Result Table:'}, ...
                'liststring', b(1, :), ...
                'SelectionMode', 'single', ...
                'ListSize', [220, 300]);
        end
        if tf ~= 1, return; end
        tab = evalin('base', a(indx).name);

    case 'Web Enrichr'
        answer = gui.myQuestdlg(FigureHandle, "Input Enrichr output table from:","Select Source", ...
            {'Paste Text', 'Open File', 'Cancel'},'Paste Text');
        switch answer
            case 'Paste Text'
            defaulttxt = sprintf('Term\tGenes\nPathway 1\tMDH1;AFMID;CAT;HYI;ACO1\nPathway 2\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nPathway 3\tGSTZ1;FAHD1;FAH;ADH5\nPathway 4\tGYS2;GBE1;PGM2\nPathway 5\tTHTPA;NFS1\n');
            % defaulttxt = sprintf('Term\tGenes\nGlyoxylate and dicarboxylate metabolism\tMDH1;AFMID;CAT;HYI;ACO1\nTryptophan metabolism\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nTyrosine metabolism\tGSTZ1;FAHD1;FAH;ADH5\nStarch and sucrose metabolism\tGYS2;GBE1;PGM2\nThiamine metabolism\tTHTPA;NFS1\nFatty acid biosynthesis\tOXSM;MCAT\nPyruvate metabolism\tMDH1;GLO1;ADH5\nPeroxisome\tSCP2;CAT;PEX1;NUDT12\nPurine metabolism\tNME7;ENTPD5;ADK;PGM2;PAICS\nPhosphonate and phosphinate metabolism\tCHPT1\nCitrate cycle (TCA cycle)\tMDH1;ACO1\nPentose phosphate pathway\tPGM2;RBKS\nbeta-Alanine metabolism\tALDH6A1;HIBCH\n');
                
            if gui.i_isuifig(FigureHandle)
                [userInput] = gui.myInputdlg({'Paste table text'}, ...
                    'Enrichr Results', {defaulttxt}, FigureHandle);
            else
                [userInput] = inputdlg('Paste table text', ...
                    'Enrichr Results', [15, 80], {defaulttxt});
            end

                if isempty(userInput)
                 disp('User canceled input.')
                 return;
                end
                a = tempname;
                fid = fopen(a, 'w');
                fprintf(fid, '%s\n', string(userInput{1}));  % Write first input only (modify for multiple)
                fclose(fid);
                tab = readtable(a,"FileType","text",'Delimiter','\t', ...
                    'VariableNamingRule', 'modify');              
            case 'Open File'
                if gui.i_isuifig(FigureHandle)
                    [fname, pathname] = uigetfile(FigureHandle, ...
                        {'*.txt', 'Enrichr Table Files (*.txt)'; ...
                        '*.*', 'All Files (*.*)'}, ...
                        'Pick an Enrichr Table File');
                else
                    [fname, pathname] = uigetfile( ...
                        {'*.txt', 'Enrichr Table Files (*.txt)'; ...
                        '*.*', 'All Files (*.*)'}, ...
                        'Pick an Enrichr Table File');
                end
                if isequal(fname, 0), return; end
                tabfile = fullfile(pathname, fname);
                warning off
                tab = readtable(tabfile, 'FileType', 'text','Delimiter','\t');                
                warning on
            otherwise
                return;
        end
        otherwise
           return;
end
end



    if all(ismember({'Term','Genes'}, tab.Properties.VariableNames))
        listitems = tab.Term;


        if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, listitems, 'Select terms:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select terms:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, 'ListSize', [260, 300]);
        end


        if tf2 ~= 1, return; end
        terms = tab.Term(indx2);
        genes = tab.Genes(indx2);
    elseif all(ismember({'TermName','OverlappingGenes'}, ...
            tab.Properties.VariableNames))
        listitems = tab.TermName;

        if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, listitems, 'Select terms:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select terms:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, 'ListSize', [260, 300]);
        end

        if tf2 ~= 1, return; end
        terms = tab.TermName(indx2);
        genes = tab.OverlappingGenes(indx2);
    else
        gui.myWarndlg(FigureHandle, 'Invalid input.');
        return;
    end
   
% taking out all gene names
allg = {};
g_by_term = {};
for ix = 1:length(genes)
    if ischar(genes{ix})
        % a = strsplit(char(genes(ix)),';'); web based
        a = strsplit(genes{ix},';');
    elseif isstring(genes{ix})
        a = strsplit(char(genes{ix}),',');
    end
    allg = [allg, a];
    g_by_term = [g_by_term, {a}]; % Get gene list for every pathway
end
allg = unique(allg);

% Creating gene counts with their corresponding pathway
df = {};
for ix = allg
    g = char(ix);
    a = {};
    for jx = 1:length(g_by_term)
        if ismember(ix, g_by_term{jx})
            a = [a,1];
        else
            a = [a,0];
        end
    end
    df = [df,a];
end

df1 = reshape(df,[length(terms),length(allg)]);
df1 = cell2table(df1,'VariableNames',allg);
df1 = table2array(df1);

hx = gui.myFigure(FigureHandle);
hx.addCustomButton('on', @in_callback_savedata, 'floppy-disk-arrow-in.jpg', 'Export data...');

CC = gui.chordChart(df1,'Arrow','off','rowName',terms,'colName',allg);
CC = CC.draw();
CC.setFont('FontName','Arial','FontSize',10)
CC.labelRotate('on');
CC.setLabelRadius(1.2);
CListC = lines(size(df1, 2));
for ix = 1:size(df1, 1)
    for jx = 1:size(df1, 2)
        CC.setChordMN(ix, jx, 'FaceColor',CListC(jx,:), 'FaceAlpha',.4)
    end
end
CC.setSquareT_Prop('FaceColor',[0,0,0])
hx.show(FigureHandle);

sz = 10;

    function in_callback_savedata(~,~)
        gui.i_exporttable(tab(indx2,:), true, ...
            'Tcircostabl','CircosTermTable',[],[], hx.FigHandle);
    end

    function in_resizefont(~, ~)
        sz = sz + 1;
        if sz > 20, sz = 5; end
        CC.setFont('FontName','Arial','FontSize', sz);
    end

    function in_PickColorMap(~, ~)
        CC.setChordColorByMap(gui.i_getrandcolormap);
    end
end