function callback_EnrichrTab2Circos(src, ~)


import mlreportgen.ppt.*;
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'resources', 'myTemplate.pptx');


    FigureHandle = src.Parent.Parent;
   
    answer = questdlg("Input Enrichr output table from:","Select Source", ...
        'Paste Text', 'Open File', 'Cancel','Paste Text');

    switch answer
        case 'Paste Text'
            [userInput] = inputdlg('Paste table text', 'Enrichr Results', ...
                [40, 20], {''});
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
            tab = readtable(tabfile);
            warning on
        otherwise
            return;
    end


    listitems = tab.Term;
    [indx2, tf2] = listdlg('PromptString', ...
    {'Select terms:'}, ...
    'SelectionMode', 'multiple', ...
    'ListString', listitems, 'ListSize', [260, 300]);

    if tf2 ~= 1, return; end

    terms = tab.Term(indx2);
    genes = tab.Genes(indx2);

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
ax0 = hFig.CurrentAxes;
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



tb = uitoolbar(hFig);
%pkg.i_addbutton2fig(tb, 'off', @in_savedata, 'export.gif', 'Export data...');
%pkg.i_addbutton2fig(tb, 'off', @in_testdata, 'plotpicker-renko.gif', 'Add Regression Line...');
%pkg.i_addbutton2fig(tb, 'off', @i_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%pkg.i_addbutton2fig(tb, 'off', @i_savemainfig2, "powerpoint.gif", 'Save Figure to PowerPoint File...');

% pkg.i_addbutton2fig(tb, 'off', @i_savemainfigx, "xpowerpoint.gif", 'Save Figure as Graphic File...');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 2}, "xpowerpoint.gif", 'Save Figure as Graphic File...');
pkg.i_addbutton2fig(tb, 'off', @i_resizefont, "icon-mat-touch-app-10.gif", 'Change Font Size');

hFig.Visible=true;

sz = 10;

    function i_resizefont(~, ~)
        sz = sz + 1;
        if sz > 20, sz = 5; end
        CC.setFont('FontName','Arial','FontSize', sz);
    end


    % function i_savemainfigx(~,~)
    %     filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
    %     [filename, filepath] = uiputfile(filter,'Save Circos Plot', ...
    %         sprintf('CircosPlot_%s','Enrichr'));
    %     if ischar(filename)
    %         exportgraphics(hFig, [filepath, filename]);
    %     end
    % end

%   function i_savemainfig2(~,~)    
%        gui.i_export2pptx({CC}, {'CircosPlot_Enrichr'});
%    end

    % function i_savemainfig(~,~)
    %     answer = questdlg('Export to PowerPoint?');
    %     if ~strcmp(answer,'Yes'), return; end
    % 
    % 
    %         OUTppt = [tempname, '.pptx'];
    %         ppt = Presentation(OUTppt, pth);
    %         open(ppt);
    %         images=cell(1,1);
    %         warning off
    %     for kx=1:1            
    %         images{kx} = [tempname, '.png'];
    %         saveas(CC,images{kx});
    %         slide3 = add(ppt, 'Small Title and Content');
    %         replace(slide3, 'Title', 'Enrichr Circos Plot');
    %         replace(slide3, 'Content', Picture(images{kx}));        
    %     end
    %         close(ppt);
    %         rptview(ppt);            
    % end

end