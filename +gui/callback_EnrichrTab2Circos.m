function callback_EnrichrTab2Circos(src, ~)

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
    'ListString', listitems, 'ListSize', [220, 300]);

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

gui.i_movegui2parent(hFig, FigureHandle);

CC=gui.chordChart(df1,'Arrow','off','rowName',terms,'colName',allg);
CC=CC.draw();

% zoom(0.75) % Zooming

% 修改字体，字号及颜色
CC.setFont('FontName','Arial','FontSize',10)
CC.labelRotate('on');

% 调节标签半径
% Adjustable Label radius
CC.setLabelRadius(1.2);    

hFig.Visible=true;
