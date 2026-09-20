function [setmatrx, setnames, setgenes, T] = enzonto(opts)
%ENZONTO  GlycoEnzOnto glycosylation pathways as a gene-set collection.
%
%   [setmatrx, setnames, setgenes] = GLY.ENZONTO() returns the
%   glycosylation pathways of GlycoEnzOnto in the same form as
%   GLY.GENESETS, so anything that takes one takes the other:
%   GLY.STATE, GLY.ENRICH, GLY.DETECT, SC_PATHWAYACTIVITY.
%
%   [setmatrx, setnames, setgenes, T] = GLY.ENZONTO() also returns a
%   table with Name, Description, Class, IsAggregate, NumGenes and Genes.
%
%   WHY, ALONGSIDE GLY.GENESETS, AND NOT INSTEAD OF IT. What this
%   buys is RESOLUTION, not coverage, and the distinction matters:
%
%     genes          GENESETS 448   GlycoEnzOnto 403
%     groupings      27 modules            95 leaf pathways (122 terms)
%     scored on one PBMC run                26 modules -> 60 pathways
%
%   So the gene list is not bigger - it is SMALLER, and differently shaped.
%   143 genes of GENESETS are absent here, and they are absent by
%   design: GlycoEnzOnto is an ontology of glycoENZYMES, so the
%   glycan-BINDING proteins are out of scope. CD22, CD207, CD209, ASGR1/2
%   and the proteoglycan core proteins ACAN, BCAN, BGN, AGRN are all in
%   GENESETS and none are here. Anything that reads a lectin module -
%   GLY.CCC's Channel A, GLY.LECTINMAP - therefore cannot be
%   switched to this collection. 98 genes go the other way, all enzymes.
%
%   Reach for this where the RESOLUTION of the prior is the binding
%   constraint, which in the scTenifoldGly work it demonstrably is: 255
%   tested ligand-receptor pairs collapsed onto 3 distinct (module, weight)
%   states with 209 sharing one. 95 pathways in place of 27 modules is a
%   direct attack on that. It is not a fix for depth confounding, which is
%   orthogonal and unmoved: on the same PBMC run the median pathway still
%   correlates 0.78 with cluster depth against 0.80 for the 27 modules.
%   Check GLY.ENRICH's rhoDepth either way.
%
%   THE FILE IS A FLATTENED ONTOLOGY, WHICH MATTERS. GlycoEnzOnto's GMT
%   lists, for every term, every gene beneath it in the hierarchy. So the
%   122 terms mix leaf pathways ("core 1 type o-N-acetylgalactosamine type
%   glycan biosynthetic pathway", 3 genes) with roll-up ancestors
%   ("glycosylation-related pathway", all 403; "glycan biosynthetic
%   pathway", 228). Scoring cells against the roll-ups is close to
%   meaningless - they saturate, and they are redundant with their own
%   children - so INCLUDE defaults to "leaf" and the 27 aggregates are held
%   back. Ask for them deliberately.
%
%   A term is called an aggregate when three or more other terms are subsets
%   of it. That is read off the gene sets rather than the OWL: it needs no
%   RDF parsing, it is exactly the property that makes a term unusable for
%   scoring, and it agrees with the ontology's own structure on this file.
%
%   CLASS is assigned the same way, by the SMALLEST aggregate that contains
%   the term: core, extension, terminal, sulfation, modification,
%   donor_synthesis, transport, degradation, regulation, or "" for the two
%   orphans (putative/unknown, and xenobiotic glucuronidation). The classes
%   reproduce the grouping Chrysinas et al. use across Tabula Sapiens (NAR
%   Genomics and Bioinformatics 2024, 6:lqae169), whose central finding is
%   that core enzymes are ubiquitous and flat across cell types while
%   terminal ones are low and cell-type specific. That makes CLASS the
%   useful axis to slice on: a module mixing core and terminal genes buries
%   the discriminative half under the invariant half.
%
%   OPTIONS:
%     opts.Include  ("leaf")  "leaf", "aggregate" or "all"
%     opts.Class    ("")      keep only these classes, e.g. ["core","terminal"]
%     opts.MinGenes (1)       drop pathways with fewer genes than this
%
%   SOURCE AND LICENCE. ASSETS/GLYCOENZONTO/GLYCOENZONTO.GMT, from
%   github.com/neel-lab/GlycoEnzOnto, CC BY 4.0, redistributed here with
%   attribution: Groth, Diehl, Gunawan and Neelamegham, "GlycoEnzOnto: a
%   GlycoEnzyme pathway and molecular function ontology", Bioinformatics 38
%   (2022) 5413-5420. The licence text is beside the file. Note the CC BY 4.0
%   is the ONTOLOGY's; the 2024 single-cell paper above is CC BY-NC and
%   nothing from it is vendored here.
%
% see also: GLY.GENESETS, GLY.STATE, GLY.ENRICH, GLY.DETECT

arguments
    opts.Include (1, 1) string {mustBeMember(opts.Include, ...
        ["leaf", "aggregate", "all"])} = "leaf"
    opts.Class (1, :) string = ""
    opts.MinGenes (1, 1) double {mustBeNonnegative} = 1
end

persistent CACHE
if isempty(CACHE)
    CACHE = i_build();
end
T = CACHE;

isagg = T.IsAggregate;
switch opts.Include
    case "leaf",      T = T(~isagg, :);
    case "aggregate", T = T(isagg, :);
end
if any(strlength(opts.Class) > 0)
    T = T(ismember(T.Class, opts.Class), :);
end
T = T(T.NumGenes >= opts.MinGenes, :);
if isempty(T)
    error("ENZONTO:Empty", ...
        "No pathway survived the Include/Class/MinGenes filters.");
end

geneCell = cell(height(T), 1);
for k = 1:height(T)
    geneCell{k} = i_split(T.Genes(k));
end
setgenes = unique(vertcat(geneCell{:}));
setgenes = setgenes(strlength(setgenes) > 0);
setnames = T.Name;
setmatrx = false(height(T), numel(setgenes));
for k = 1:height(T)
    setmatrx(k, :) = ismember(upper(setgenes), upper(geneCell{k}));
end

end


% ----------------------------------------------------------------------
function T = i_build()
here = fileparts(fileparts(mfilename('fullpath')));
f = fullfile(here, 'assets', 'GlycoEnzOnto', 'GlycoEnzOnto.gmt');
if ~isfile(f)
    error("ENZONTO:NoAsset", ...
        ['GlycoEnzOnto.gmt not found at %s. Run ', ...
        'ASSETS/GLYCOENZONTO/FETCH_GLYCOENZONTO to download it.'], f);
end

lines = readlines(f);
lines = lines(strlength(strtrim(lines)) > 0);
n = numel(lines);
name = strings(n, 1); desc = strings(n, 1); genes = cell(n, 1);
for k = 1:n
    p = split(lines(k), sprintf('\t'));
    name(k) = i_unquote(p(1));
    if numel(p) > 1, desc(k) = i_unquote(p(2)); end
    if numel(p) > 2
        g = arrayfun(@i_unquote, p(3:end));
        genes{k} = unique(upper(g(strlength(g) > 0)));
    else
        genes{k} = strings(0, 1);
    end
end

% Containment, once, over all pairs. A term with three or more other terms
% wholly inside it is a roll-up rather than a pathway to score against.
contains_n = zeros(n, 1);
for i = 1:n
    if isempty(genes{i}), continue, end
    for j = 1:n
        if i == j || isempty(genes{j}), continue, end
        if all(ismember(genes{j}, genes{i})) && ...
                numel(genes{j}) < numel(genes{i})
            contains_n(i) = contains_n(i) + 1;
        end
    end
end
isAgg = contains_n >= 3;

% Class by the smallest aggregate that contains the term.
classOf = i_classroots();
cls = strings(n, 1);
rootIdx = zeros(numel(classOf), 1);
for r = 1:numel(classOf)
    hit = find(name == classOf(r).name, 1);
    if ~isempty(hit), rootIdx(r) = hit; end
end
for k = 1:n
    best = ""; bestSize = Inf;
    for r = 1:numel(classOf)
        ri = rootIdx(r);
        % A term that IS a class root carries that class itself; the
        % smallest-containing rule then still resolves its children.
        if ri == 0 || isempty(genes{k}), continue, end
        if all(ismember(genes{k}, genes{ri})) && numel(genes{ri}) <= bestSize
            best = classOf(r).class;
            bestSize = numel(genes{ri});
        end
    end
    cls(k) = best;
end

joined = strings(n, 1);
for k = 1:n
    joined(k) = strjoin(genes{k}, ",");
end
T = table(name, desc, cls, isAgg, cellfun(@numel, genes), joined, ...
    VariableNames = ["Name", "Description", "Class", "IsAggregate", ...
    "NumGenes", "Genes"]);
T = sortrows(T, ["IsAggregate", "Class", "Name"]);
end


function c = i_classroots()
% The ontology's own mid-level terms, used as class labels. Ordered here for
% readability only; assignment picks the smallest containing one.
spec = { ...
    "glycan core structure biosynthetic pathway",           "core"; ...
    "glycan elongation and branching biosynthetic pathway", "extension"; ...
    "glycan capping structure biosynthetic pathway",        "terminal"; ...
    "glycan monosaccharide sulfation pathway",              "sulfation"; ...
    "glycan monosaccharide modification pathway",           "modification"; ...
    "glycan molecular donor synthesis pathway",             "donor_synthesis"; ...
    "glycan substrate transport pathway",                   "transport"; ...
    "glycan degradation pathway",                           "degradation"; ...
    "glycan biosynthesis regulation",                       "regulation"};
c = struct('name', spec(:, 1), 'class', spec(:, 2));
end


function s = i_unquote(s)
s = strtrim(string(s));
if strlength(s) >= 2 && startsWith(s, '"') && endsWith(s, '"')
    s = extractBetween(s, 2, strlength(s) - 1);
end
s = string(s);
end


function v = i_split(s)
v = strtrim(split(string(s), ","));
v = v(strlength(v) > 0);
end
