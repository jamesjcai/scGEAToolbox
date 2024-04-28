function callback_EnrichrTab2Circos(src, ~)


import mlreportgen.ppt.*;
%pw1 = fileparts(mfilename('fullpath'));
%pth = fullfile(pw1, '..', 'resources', 'myTemplate.pptx');


    FigureHandle = src.Parent.Parent;
   
    answer = questdlg("Input Enrichr output table from:","Select Source", ...
        'Paste Text', 'Open File', 'Cancel','Paste Text');

    switch answer
        case 'Paste Text'
defaulttxt = sprintf('Term\tGenes\nPathway a\tMDH1;AFMID;CAT;HYI;ACO1\nPathway b\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nPathway c\tGSTZ1;FAHD1;FAH;ADH5\nPathway d\tGYS2;GBE1;PGM2\nPathway e\tTHTPA;NFS1\n');
% defaulttxt = sprintf('Term\tGenes\nGlyoxylate and dicarboxylate metabolism\tMDH1;AFMID;CAT;HYI;ACO1\nTryptophan metabolism\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nTyrosine metabolism\tGSTZ1;FAHD1;FAH;ADH5\nStarch and sucrose metabolism\tGYS2;GBE1;PGM2\nThiamine metabolism\tTHTPA;NFS1\nFatty acid biosynthesis\tOXSM;MCAT\nPyruvate metabolism\tMDH1;GLO1;ADH5\nPeroxisome\tSCP2;CAT;PEX1;NUDT12\nPurine metabolism\tNME7;ENTPD5;ADK;PGM2;PAICS\nPhosphonate and phosphinate metabolism\tCHPT1\nCitrate cycle (TCA cycle)\tMDH1;ACO1\nPentose phosphate pathway\tPGM2;RBKS\nbeta-Alanine metabolism\tALDH6A1;HIBCH\n');
            
            [userInput] = inputdlg('Paste table text', 'Enrichr Results', ...
                [15, 60], {defaulttxt});
            if isempty(userInput)
             disp('User canceled input.')
             return;
            end
            a = tempname;
            fid = fopen(a, 'w');
            fprintf(fid, '%s\n', string(userInput{1}));  % Write first input only (modify for multiple)
            fclose(fid);
            warning off
            tab = readtable(a,"FileType","text");
            warning on
        case 'Open File'
            [fname, pathname] = uigetfile( ...
                {'*.txt', 'Enrichr Table Files (*.txt)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick an Enrichr Table File');
            if isequal(fname, 0), return; end
            tabfile = fullfile(pathname, fname);
            warning off
            tab = readtable(tabfile, 'FileType', 'text');
            warning on
        otherwise
            return;
    end


    if all(ismember({'Term','Genes'}, tab.Properties.VariableNames))
        listitems = tab.Term;    
        [indx2, tf2] = listdlg('PromptString', ...
        {'Select terms:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, 'ListSize', [260, 300]);
        if tf2 ~= 1, return; end
        terms = tab.Term(indx2);
        genes = tab.Genes(indx2);
    else
        listitems = tab.(tab.Properties.VariableNames{1});
        [indx2, tf2] = listdlg('PromptString', ...
        {'Select terms:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, 'ListSize', [260, 300]);
        if tf2 ~= 1, return; end
        terms = tab.(tab.Properties.VariableNames{1})(indx2);
        genes = tab.(tab.Properties.VariableNames{2})(indx2);
    end

% taking out all gene names
allg = {};
g_by_term = {};
for i = 1:length(genes)
    a = strsplit(char(genes(i)),';');
    allg = [allg,a];
    g_by_term = [g_by_term,{a}]; % Get gene list for every pathway
end
allg = unique(allg);

% Creating gene counts with their corresponding pathway
df = {};
for i = allg
    g = char(i);
    a = {};
    for j = 1:length(g_by_term)
        if ismember(i,g_by_term{j})
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

hFig = figure('Visible',"off");
% ax0 = hFig.CurrentAxes;
gui.i_movegui2parent(hFig, FigureHandle);

CC = gui.chordChart(df1,'Arrow','off','rowName',terms,'colName',allg);
CC = CC.draw();
% zoom(0.75) % Zooming
% 修改字体，字号及颜色
CC.setFont('FontName','Arial','FontSize',10)
CC.labelRotate('on');
% 调节标签半径
% Adjustable Label radius
CC.setLabelRadius(1.2);

% 修改弦颜色(Modify chord color)
CListC = lines(size(df1, 2));
for i = 1:size(df1, 1)
    for j = 1:size(df1, 2)
        CC.setChordMN(i,j, 'FaceColor',CListC(j,:), 'FaceAlpha',.4)
    end
end

tb = findall(hFig, 'Tag', 'FigureToolBar');
uipushtool(tb, 'Separator', 'off');

% tb = uitoolbar(hFig);
%pkg.i_addbutton2fig(tb, 'off', @in_savedata, 'export.gif', 'Export data...');
%pkg.i_addbutton2fig(tb, 'off', @in_testdata, 'plotpicker-renko.gif', 'Add Regression Line...');
%pkg.i_addbutton2fig(tb, 'off', @i_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%pkg.i_addbutton2fig(tb, 'off', @i_savemainfig2, "powerpoint.gif", 'Save Figure to PowerPoint File...');

% pkg.i_addbutton2fig(tb, 'off', @i_savemainfigx, "xpowerpoint.gif", 'Save Figure as Graphic File...');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 2}, "xpowerpoint.gif", 'Save Figure as Graphic File...');
pkg.i_addbutton2fig(tb, 'off', @i_resizefont, "noun_font_size_591141.gif", 'Change Font Size');
pkg.i_addbutton2fig(tb, 'off', @i_PickColorMap, "plotpicker-compass.gif", "Pick new color map");


hFig.Visible=true;
sz = 10;

    function i_resizefont(~, ~)
        sz = sz + 1;
        if sz > 20, sz = 5; end
        CC.setFont('FontName','Arial','FontSize', sz);
    end

    function i_PickColorMap(~, ~)
        CC.setChordColorByMap(gui.i_getrandcolormap);
    end

end