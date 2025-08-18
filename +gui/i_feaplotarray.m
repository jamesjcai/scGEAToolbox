function i_feaplotarray(sce, tgene, thisc, uselog, parentfig)

if nargin < 5, parentfig = []; end
if nargin < 4, uselog = false; end

[Xt] = gui.i_transformx(sce.X, [], [], parentfig);
if isempty(Xt), return; end

X = Xt;
g = sce.g;
s = sce.s;

if max(findgroups(thisc))>50
    answer = gui.myQuestdlg(parentfig, 'Too many groups. Continue?','',[],[],'warning');
    if ~strcmp(answer, 'Yes'), return; end
end

[c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
if noanswer, return; end

cL = strrep(cL(:),'_','\_');
[yes] = ismember(tgene, g);
if ~any(yes), warning('No genes found.'); return; end
z = length(tgene) - sum(yes);
if z > 0, fprintf('%d gene(s) not in the list are excluded.\n', z); end
tgene = tgene(yes);
if issparse(X), X = full(X); end
if uselog, X = log1p(X); end

a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

%z = zeros(size(X, 2)*length(tgene),1);
z = [];
for kx = 1:length(tgene)
    z =[z, X(g == tgene(kx), :)];
end

for kx = 1:length(tgene)

    hx=gui.myFigure;
    %hx = figure;
    for ky = 1:length(cL)
        cellidx = c==ky;
        nexttile;

        % sc_scattermarker(Xt, sce.g, sce.s, glist(k), 1, 5, false);   

        ydata = X(g == tgene(kx), cellidx);
    
        if size(s,2)>2
            scatter3(s(cellidx, 1), s(cellidx, 2), s(cellidx, 3), 5, ydata, 'filled');
        else
            scatter(s(cellidx, 1), s(cellidx, 2), 5, ydata, 'filled');
        end
        gui.i_setautumncolor(ydata, a, true, any(ydata==0));
        clim([min(z) max(z)]);  % Adjust color axis to data range
        title(cL{ky});
    end
    sgtitle(tgene(kx));
    hx.show(parentfig);
end
