function [Tm, Tw] = i_markerlist2weight(sce, FigureHandle)
if nargin<2, FigureHandle = []; end
%Tm=readtable('markerlist_hs.txt','ReadVariableNames',false);
%celltypev=string(Tm.Var1);
%markergenev=string(Tm.Var2);
Tm = [];
Tw = [];
if nargin < 1, sce = []; end

answer = gui.myQuestdlg(FigureHandle, 'Load scType marker gene list?');
indata = '';
if strcmp(answer, 'Yes')
    indata = gui.i_getsctypemarkers(FigureHandle);
end

if isempty(indata)
    if isempty(sce)
        indata = sprintf('Cell type 1\tgene1,gene2,gene3,gene4\nCell type 2\tgene5,gene6,gene7');
    else
        % try
        %     t=sc_splinefit(sce.X,sce.g);
        %     a=t.genes;
        % catch ME
        %     warning(ME.message);
        %     a = sce.g(randperm(length(sce.g)));
        % end
        rng("shuffle");
        a = sce.g(randperm(length(sce.g)));
        a1 = sprintf('%s,%s,%s,%s', a(1), a(2), a(3), a(4));
        a2 = sprintf('%s,%s,%s', a(5), a(6), a(7));
        indata = sprintf('Cell type 1\t%s\nCell type 2\t%s', a1, a2);
    end
end

% indata=gui.i_getsctypemarkers;
% assignin("base","indata",indata);

if gui.i_isuifig(FigureHandle)

    a = gui.myInputwin([], [], indata, FigureHandle);
    % a = gui.myInputdlg({sprintf('Format:\nCell type name [TAB] Gene1,Gene2')}, ...
    %     'Markers Input', {char(indata)}, FigureHandle);
 else
    a = inputdlg(sprintf('Format:\nCell type name [TAB] Gene1,Gene2'), ...
        'Markers Input', [15, 80], {char(indata)}, 'on');
 end


if isempty(a), return; end
b = strtrim(string(a{1}));
[c, d] = strtok(b, sprintf('\t'));
d = upper(d);
d = strrep(d, ' ', '');
Tm = table(strtrim(c), strtrim(d));


s = upper(string(Tm.Var2));
S = [];
for k = 1:length(s)
    a = strsplit(s(k), ',');
    a = strtrim(a);
    if strlength(a(end)) == 0 || isempty(a(end))
        a = a(1:end-1);
    end
    S = [S, a];
end

%%
N = length(S);
t = tabulate(S);
f = cell2mat(t(:, 3));
if max(f) - min(f) < eps
    w = ones(N, 1);
else
    w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
end
genelist = string(t(:, 1));
Tw = table(genelist, w);
Tw.Properties.VariableNames = {'Var1', 'Var2'};
