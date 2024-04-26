% Enrichr Circos Plots

% Read in Enrichr Table
tab = readtable('KEGG_2021_Human_table.txt');

% Select terms
Termlist = {'Glyoxylate and dicarboxylate metabolism','Tryptophan metabolism'};

% Subset the table using desired terms
% tab_sub = tab(1:5,:);
tab_sub = tab(ismember(tab.Term,Termlist),:);

% Extracting terms
terms = tab_sub.Term;

% Extracting genes
genes = tab_sub.Genes;

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

figure();
CC=gui.chordChart(df1,'Arrow','off','rowName',terms,'colName',allg);
CC=CC.draw();

zoom(0.75) % Zooming

% 修改字体，字号及颜色
CC.setFont('FontName','Arial','FontSize',10)
CC.labelRotate('on');

% 调节标签半径
% Adjustable Label radius
CC.setLabelRadius(1.2);