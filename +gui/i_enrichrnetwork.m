function i_enrichrnetwork()

defaulttxt = sprintf('Term\tGenes\nPathway 1\tMDH1;AFMID;CAT;HYI;ACO1\nPathway 2\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nPathway 3\tGSTZ1;FAHD1;FAH;ADH5\nPathway 4\tGYS2;GBE1;PGM2\nPathway 5\tTHTPA;NFS1\n');
% defaulttxt = sprintf('Term\tGenes\nGlyoxylate and dicarboxylate metabolism\tMDH1;AFMID;CAT;HYI;ACO1\nTryptophan metabolism\tAFMID;CAT;KMO;ALDH8A1;DHTKD1\nTyrosine metabolism\tGSTZ1;FAHD1;FAH;ADH5\nStarch and sucrose metabolism\tGYS2;GBE1;PGM2\nThiamine metabolism\tTHTPA;NFS1\nFatty acid biosynthesis\tOXSM;MCAT\nPyruvate metabolism\tMDH1;GLO1;ADH5\nPeroxisome\tSCP2;CAT;PEX1;NUDT12\nPurine metabolism\tNME7;ENTPD5;ADK;PGM2;PAICS\nPhosphonate and phosphinate metabolism\tCHPT1\nCitrate cycle (TCA cycle)\tMDH1;ACO1\nPentose phosphate pathway\tPGM2;RBKS\nbeta-Alanine metabolism\tALDH6A1;HIBCH\n');

[userInput] = inputdlg('Paste table text', 'Enrichr Results', ...
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
tab = readtable(a,"FileType","text",'Delimiter','\t');                
warning on

terms = tab.Term;
genes = tab.Genes;

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
df1 = reshape(df, [length(terms), length(allg)]);
df1 = cell2table(df1, 'VariableNames', allg);
df1 = table2array(df1);

A = zeros(length(terms)+length(allg));
A(1:height(df1), length(terms)+1:end) = df1;
A = A+A';
names = [terms; allg'];
G = graph(A, names);
% layoutcoords(G,"force");

net = gui.networkvisualizer(A);
net.setNodeLabels(names);
net.setNodeSizes('auto');
% figure; plot(net);
%%
fig = figure;
p = plot(G, 'ButtonDownFcn', @startDragFcn);
p.layout("force");

ix = 1:length(terms);
if any(ix)
    cc = repmat([0, 0, 0], G.numnodes, 1);
    cc(ix, :) = repmat([1, 0, 0], length(terms), 1);
    p.NodeLabelColor = cc;
end

    % Callback to initiate dragging
    function startDragFcn(hObj, ~)
        
        %fig = ancestor(hObj, 'figure');
        set(fig, 'WindowButtonMotionFcn', {@draggingFcn, hObj});
        set(fig, 'WindowButtonUpFcn', @stopDragFcn);
    end
    
    % Function to drag the point
    function draggingFcn(~, ~, hObj)
        % Current cursor position in data coordinates
        cp = get(gca, 'CurrentPoint');
        % Update the y-data of the nearest point
        yData = get(hObj, 'YData');
        xData = get(hObj, 'XData');
        % xData
        
        idx = dsearchn([xData' yData'], [cp(1,1) cp(1,2)]);
        %[~, idx] = min(abs(xData - cp(1,1)));   % Find closest x to mouse
        xData(idx) = cp(1,1);
        yData(idx) = cp(1,2);     
        set(hObj, 'XData', xData); % Update y value
        set(hObj, 'YData', yData);
    end
    
    % Function to stop dragging
    function stopDragFcn(~, ~)
        %fig = gcbf;
        set(fig, 'WindowButtonMotionFcn', '');
        set(fig, 'WindowButtonUpFcn', '');
    end

    end