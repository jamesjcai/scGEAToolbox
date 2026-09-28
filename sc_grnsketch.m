function [hFig, info] = sc_grnsketch(A, genes, options)
% SC_GRNSKETCH  Sketch a genome-wide GRN as a map of gene modules.
%
%   sc_grnsketch(A, genes) draws the network A (genes x genes, e.g. from
%   net.pcrnet) as a map of its gene modules. A few thousand genes cannot
%   be drawn legibly as genes, so the sketch has two levels:
%
%     Module map  one node per module, sized by its gene count and labelled
%                 with its top hub genes; an edge's width is the mean
%                 |weight| between the two modules. Modules are Louvain
%                 communities of the graph that keeps each gene's
%                 Neighbors strongest partners.
%     Gene view   click a module to open its NumNodes strongest hub genes
%                 with their strongest edges, or use the Find Gene button
%                 to open one gene with its strongest partners.
%
%   sc_grnsketch(file) reads A and the gene list (g, genes or genelist)
%   from a .mat file, e.g. one saved by Network > Build Gene Regulatory Network (GRN)... (All Genes).
%
%   [hFig, info] = sc_grnsketch(...) also returns
%     info.Module       Nx1 module index per gene, 1 = largest; 0 for genes
%                       with no edges
%     info.Hub          Nx1 within-module strength (sum of |weight| to the
%                       gene's own module), the hub ranking
%     info.Modules      table: Module, NumGenes, Cohesion (mean within-module
%                       |weight|), TopHubs
%     info.ModuleWeight MxM mean |weight| between modules
%
%   Name-value arguments
%     NumModules    Modules drawn on the map, largest first (default 20)
%     MinModuleSize Smaller modules are not drawn (default 10)
%     Neighbors     Strongest partners kept per gene for module detection
%                   (default 10)
%     Resolution    Louvain resolution; larger gives more, smaller modules.
%                   Default: resolution 1, then the largest modules are split
%                   one at a time, until min(12, NumModules) modules reach
%                   MinModuleSize. A split is kept only when its parts are
%                   at least 1.25x more tightly connected within than
%                   between, so a module with no internal structure stays
%                   whole. info.Modules.Properties.UserData records the
%                   resolution and the number of splits.
%     NumNodes      Genes drawn in a gene view (default 60)
%     EdgeCutoff    Quantile of |weight| below which gene-view edges are
%                   hidden (default 0.9)
%     MinDegree     Every gene-view node keeps its MinDegree strongest
%                   edges (default 2)
%     Module        Open this module's gene view directly, without the map
%     FocusGene     Open this gene's neighbourhood directly, without the map
%     ParentFig     Parent figure for theme and placement
%
%   Example
%     [A, g] = deal(net.pcrnet(X), genelist);
%     [~, info] = sc_grnsketch(A, g);
%     sc_grnsketch(A, g, FocusGene="STAT1");
%
% See also sc_grnview, ten.sctenifoldnetview, pkg.e_louvain, net.pcrnet.

arguments
    A
    genes = []
    options.NumModules (1, 1) double {mustBePositive, mustBeInteger} = 20
    options.MinModuleSize (1, 1) double {mustBePositive, mustBeInteger} = 10
    options.Neighbors (1, 1) double {mustBePositive, mustBeInteger} = 10
    options.Resolution double {mustBeScalarOrEmpty, mustBePositive} = []
    options.NumNodes (1, 1) double {mustBePositive, mustBeInteger} = 60
    options.EdgeCutoff (1, 1) double {mustBeGreaterThanOrEqual(options.EdgeCutoff, 0), ...
        mustBeLessThan(options.EdgeCutoff, 1)} = 0.9
    options.MinDegree (1, 1) double {mustBeNonnegative, mustBeInteger} = 2
    options.Module double = []
    options.FocusGene string = string.empty
    options.ParentFig = []
end

[A, genes] = i_resolveinput(A, genes);
n = size(A, 1);
A = double(full(A));
A = 0.5*(A + A.');
A(1:(n + 1):end) = 0;
absA = abs(A);

if ~isempty(options.FocusGene)
    hFig = i_geneview(A, genes, options.FocusGene, options, options.ParentFig);
    info = struct;
    return;
end

[module, hub, T, W] = i_modules(absA, genes, options);
info = struct("Module", module, "Hub", hub, "Modules", T, "ModuleWeight", W);

if ~isempty(options.Module)
    if ~isscalar(options.Module) || ~ismember(options.Module, T.Module)
        error("sc_grnsketch:badModule", ...
            "Module must be one of 1..%d.", height(T));
    end
    hFig = i_moduleview(A, genes, module, hub, options.Module, options, options.ParentFig);
    return;
end

hFig = i_modulemap(A, genes, module, hub, T, W, options);
end

% ------------------------------------------------------------------ input

function [A, genes] = i_resolveinput(A, genes)
if ischar(A) || (isstring(A) && isscalar(A))
    fname = char(A);
    if ~isfile(fname)
        error("sc_grnsketch:fileNotFound", "File not found: %s", fname);
    end
    S = load(fname);
    if ~isfield(S, "A")
        error("sc_grnsketch:noAdjacency", "%s has no variable named A.", fname);
    end
    A = S.A;
    if isempty(genes)
        for name = ["g", "genes", "genelist"]
            if isfield(S, name)
                genes = S.(name);
                break
            end
        end
    end
end
if ~(isnumeric(A) || islogical(A)) || ~ismatrix(A) || size(A, 1) ~= size(A, 2)
    error("sc_grnsketch:badAdjacency", "A must be a square adjacency matrix.");
end
if isempty(genes)
    genes = "G" + string(1:size(A, 1)).';
end
genes = string(genes(:));
if numel(genes) ~= size(A, 1)
    error("sc_grnsketch:geneCount", ...
        "genes has %d entries but A is %dx%d.", numel(genes), size(A, 1), size(A, 2));
end
end

% ---------------------------------------------------------------- modules

function [module, hub, T, W] = i_modules(absA, genes, options)
% Louvain on the union of every gene's strongest partners. A dense PCR
% network has an edge between every pair, so modularity on it directly is
% dominated by weak background weights; keeping each gene's strongest
% partners is the usual kNN sparsification and keeps the graph connected.
n = size(absA, 1);
active = any(absA > 0, 2);
k = min(options.Neighbors, n - 1);
[val, col] = maxk(absA, k, 2);
keep = val > 0;
rows = repmat((1:n).', 1, k);
K = sparse(rows(keep), col(keep), val(keep), n, n);
K = max(K, K.');

raw = zeros(n, 1);
[raw(active), gamma, nSplit] = i_louvain(K(active, active), absA(active, active), options);

% Relabel by size, largest first; 0 stays "no edges"
labels = unique(raw(raw > 0));
sizes = accumarray(raw(raw > 0), 1);
sizes = sizes(labels);
[~, order] = sort(sizes, "descend");
module = zeros(n, 1);
for m = 1:numel(order)
    module(raw == labels(order(m))) = m;
end
nModule = numel(order);

% Strength of every gene towards every module, from the full network
member = sparse(find(module > 0), module(module > 0), 1, n, nModule);
toModule = absA*member;
hub = zeros(n, 1);
hub(module > 0) = toModule(sub2ind(size(toModule), find(module > 0), module(module > 0)));

sizes = full(sum(member, 1)).';
total = full(member.'*toModule);
pairs = sizes*sizes.';
pairs(1:(nModule + 1):end) = sizes.*(sizes - 1);
W = total./max(pairs, 1);
cohesion = diag(W);
W(1:(nModule + 1):end) = 0;

nTop = 3;
topHubs = strings(nModule, 1);
for m = 1:nModule
    idx = find(module == m);
    [~, o] = sort(hub(idx), "descend");
    topHubs(m) = strjoin(genes(idx(o(1:min(nTop, numel(o))))), ", ");
end
T = table((1:nModule).', sizes, cohesion, topHubs, ...
    VariableNames=["Module", "NumGenes", "Cohesion", "TopHubs"]);
T.Properties.UserData = struct("Resolution", gamma, "NumSplits", nSplit);
end

function [c, gamma, nSplit] = i_louvain(K, absA, options)
% Louvain at the given resolution, or, by default, at resolution 1 followed
% by supported splits of the largest modules.
%
% A PCR network is built from a few components, so at resolution 1 its
% genes split little further than the components they load on, which is
% too coarse to sketch. Raising the resolution globally splits everything,
% including modules with no internal structure. Instead each large module
% is split on its own, and the split is kept only if its parts are clearly
% more tightly connected within than between on the full network - not on
% the kNN graph, where selecting the strongest edges makes any partition
% look good.
minContrast = 1.25;
nSplit = 0;
if ~isempty(options.Resolution)
    gamma = options.Resolution;
    c = pkg.e_louvain(K, gamma);
    return;
end
gamma = 1;
c = pkg.e_louvain(K, gamma);
target = min(12, options.NumModules);
minSize = options.MinModuleSize;
unsplittable = false(max(c), 1);
while true
    sizes = accumarray(c(:), 1);
    if nnz(sizes >= minSize) >= target, break; end
    sizes(unsplittable | sizes < 2*minSize) = 0;
    [s, m] = max(sizes);
    if s == 0, break; end
    idx = find(c == m);
    sub = pkg.e_louvain(K(idx, idx), 1);
    if max(sub) > 1 && i_splitcontrast(absA(idx, idx), sub) >= minContrast
        c(idx(sub > 1)) = max(c) + sub(sub > 1) - 1;
        unsplittable(end + 1:max(c)) = false;
        nSplit = nSplit + 1;
    else
        unsplittable(m) = true;
    end
end
end

function r = i_splitcontrast(absA, sub)
% Mean |weight| within the parts over mean |weight| between them
same = sub(:) == sub(:).';
offdiag = ~eye(numel(sub));
within = mean(absA(same & offdiag));
between = mean(absA(~same));
r = within/max(between, realmin);
end

% ------------------------------------------------------------- module map

function hFig = i_modulemap(A, genes, module, hub, T, W, options)
drawn = find(T.NumGenes >= options.MinModuleSize, options.NumModules);
if isempty(drawn)
    drawn = 1:min(options.NumModules, height(T));
end
nd = numel(drawn);
Wd = W(drawn, drawn);

% Keep the strongest quarter of module links, and each module's strongest
% one, so no module floats free without saying where it attaches.
Wd = i_prune(Wd, i_quantile(Wd, 0.75), 1);
names = "M" + T.Module(drawn) + " (" + T.NumGenes(drawn) + ")";
G = graph(Wd, cellstr(names), "omitselfloops");

hx = gui.myFigure(options.ParentFig);
hFig = hx.FigHandle;
ax = hx.AxHandle;
if G.numedges > 0
    lw = 0.5 + 4.5*G.Edges.Weight/max(G.Edges.Weight);
else
    lw = 0.5;
end
p = plot(ax, G, "Layout", "force", "WeightEffect", "inverse", ...
    "LineWidth", lw, "EdgeColor", [0.6, 0.6, 0.6], "EdgeAlpha", 0.6);
colors = i_modulecolors(nd);
p.NodeColor = colors;
p.MarkerSize = 6 + 24*sqrt(T.NumGenes(drawn)/max(T.NumGenes(drawn)));
p.NodeLabel = cellstr(names + ": " + T.TopHubs(drawn));
p.NodeFontSize = 8;
p.Interpreter = "none";
p.ButtonDownFcn = @in_click;
% Hide only the rulers, keeping the themed axes background, and paint the
% figure white in light mode so the map does not sit on theme grey.
ax.XAxis.Visible = "off";
ax.YAxis.Visible = "off";
box(ax, "off");
i_whitebackground(hFig);
if isprop(hFig, "ThemeChangedFcn")
    hFig.ThemeChangedFcn = @(src, ~) i_whitebackground(src);
end
nHidden = height(T) - nd;
title(ax, sprintf("%d genes in %d modules (Louvain resolution %g, %d supported splits)", ...
    nnz(module), height(T), T.Properties.UserData.Resolution, ...
    T.Properties.UserData.NumSplits), "Interpreter", "none");
subtitle(ax, sprintf(['Node size: genes in the module (top hub genes listed). ', ...
    'Edge width: mean |weight| between modules.\nClick a module to open its ', ...
    '%d strongest hub genes.%s'], options.NumNodes, ...
    i_hiddennote(nHidden, options.MinModuleSize)), "Interpreter", "none");

hx.addCustomButton("off", @in_findgene, "network_node_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", ...
    "Find Gene: open one gene with its strongest partners");
hx.addCustomButton("off", @in_table, "data_table_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", ...
    "Show Module Table");
hx.addCustomButton("off", @in_export, "export.gif", "Export Modules to Workspace");
hx.setTitle("GRN Sketch");
hx.show(options.ParentFig);

    function in_click(~, ~)
        cp = get(ax, "CurrentPoint");
        k = dsearchn([p.XData(:), p.YData(:)], cp(1, 1:2));
        gui.myFigure.drawInto(hFig, ...
            @() i_moduleview(A, genes, module, hub, drawn(k), options, hFig));
    end

    function in_findgene(~, ~)
        answer = gui.myInputdlg({'Gene name:'}, 'Find Gene', {char(extractBefore(T.TopHubs(1) + ",", ","))}, hFig);
        if isempty(answer) || isempty(strtrim(answer{1})), return; end
        try
            gene = string(strtrim(answer{1}));
            gui.myFigure.drawInto(hFig, @() i_geneview(A, genes, gene, options, hFig));
        catch ME
            gui.myErrordlg(hFig, ME.message, ME.identifier);
        end
    end

    function in_table(~, ~)
        % In the sketch's place on screen, with a Back button to it.
        gui.i_openinstead(hFig, @() gui.TableViewerApp(T, hFig));
    end

    function in_export(~, ~)
        if ~(ismcc || isdeployed)
            export2wsdlg({'Save module table to variable named:', ...
                'Save module index per gene to variable named:'}, ...
                {'grnmodules', 'genemodule'}, {T, module});
        end
    end
end

function i_whitebackground(hFig)
% White in light mode. A manually set Color stays white when the theme
% turns dark, and neither resetting ColorMode nor re-applying the same
% theme repaints it, so the theme's own colour is read off a hidden figure.
style = "light";
try
    style = string(hFig.Theme.BaseColorStyle);
catch
    % no Theme before R2025a; the classic figure is light
end
if style == "light"
    hFig.Color = "w";
else
    tmp = figure(Visible="off");
    theme(tmp, style);
    hFig.Color = tmp.Color;
    delete(tmp);
end
end

function s = i_hiddennote(nHidden, minSize)
if nHidden > 0
    s = sprintf(" %d smaller modules (under %d genes or past the limit) are not drawn.", ...
        nHidden, minSize);
else
    s = "";
end
end

% -------------------------------------------------------------- gene views

function hFig = i_moduleview(A, genes, module, hub, m, options, parentfig)
idx = find(module == m);
[~, o] = sort(hub(idx), "descend");
idx = idx(o(1:min(options.NumNodes, numel(o))));
nHub = min(5, numel(idx));
highlight = false(numel(idx), 1);
highlight(1:nHub) = true;
figname = sprintf("Module M%d: top %d of %d genes by hub strength", ...
    m, numel(idx), nnz(module == m));
hFig = i_drawgenes(A, genes, idx, highlight, sprintf("top %d hub genes", nHub), ...
    "module genes", figname, options, parentfig);
end

function hFig = i_geneview(A, genes, gene, options, parentfig)
g0 = find(strcmpi(genes, gene), 1);
if isempty(g0)
    error("sc_grnsketch:geneNotFound", "%s is not in the network.", gene);
end
[val, o] = sort(abs(A(g0, :)), "descend");
o = o(val > 0);
idx = [g0; o(1:min(options.NumNodes - 1, numel(o))).'];
highlight = false(numel(idx), 1);
highlight(1) = true;
figname = sprintf("%s and its %d strongest partners", genes(g0), numel(idx) - 1);
hFig = i_drawgenes(A, genes, idx, highlight, "focus gene", "partners", ...
    figname, options, parentfig);
end

function hFig = i_drawgenes(A, genes, idx, highlight, hname, pname, figname, options, parentfig)
sub = A(idx, idx);
t = i_quantile(sub, options.EdgeCutoff);
sub = i_prune(sub, t, options.MinDegree);
G = graph(sub, cellstr(genes(idx)), "omitselfloops");
nodeinfo = struct;
nodeinfo.PanelTitles = string(sprintf("%d genes, %d edges", numel(idx), G.numedges));
nodeinfo.Highlight = highlight;
nodeinfo.HighlightName = hname;
nodeinfo.PlainName = pname;
nodeinfo.IsTF = i_istf(genes(idx));
nodeinfo.EdgeRef = max(abs(sub(:)));
nodeinfo.Legend = sprintf(['\nEdges below |weight| %.3g (the %g quantile among ', ...
    'these genes) are hidden, but every gene keeps its %d strongest edges.'], ...
    t, options.EdgeCutoff, options.MinDegree);
hFig = gui.i_multigraphs({G}, nodeinfo, char(figname), parentfig);
end

% ---------------------------------------------------------------- helpers

function t = i_quantile(M, q)
a = abs(nonzeros(triu(M, 1)));
if isempty(a)
    t = 0;
else
    t = quantile(a, q);
end
end

function M = i_prune(M, t, mindeg)
% Hide edges under t, but let each node keep its mindeg strongest edges,
% so the sketch never shows a node without saying where it attaches.
M = M - diag(diag(M));
M = 0.5*(M + M.');
keep = abs(M) >= t & M ~= 0;
if mindeg > 0
    [~, order] = sort(abs(M), 2, "descend");
    for k = 1:size(M, 1)
        top = order(k, 1:min(mindeg, size(M, 2)));
        top = top(M(k, top) ~= 0);
        keep(k, top) = true;
        keep(top, k) = true;
    end
end
M = M.*keep;
end

function c = i_modulecolors(n)
base = orderedcolors("gem12");
c = base(mod((0:n - 1).', size(base, 1)) + 1, :);
end

function tf = i_istf(genelist)
tf = false(numel(genelist), 1);
fname = fullfile(fileparts(mfilename("fullpath")), "assets", "TFome", "tfome_tfgenes.mat");
if ~isfile(fname), return; end
S = load(fname, "tfgenes");
tf = ismember(upper(string(genelist)), string(S.tfgenes));
end
