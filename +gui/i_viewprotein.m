function info = i_viewprotein(query, options)
%I_VIEWPROTEIN Show the 3-D structure of a gene's protein, a PDB entry, or a file.
%   INFO = gui.i_viewprotein(QUERY) opens a structure viewer for QUERY,
%   which is one of:
%     - a gene symbol, e.g. "EGFR" or "Egfr"
%     - a four-character PDB ID, e.g. "1IYT"
%     - a local structure file (.pdb, .pdbqt, .cif, .sdf, .mol2, ...)
%
%   A gene is resolved in this order: the curated gene -> PDB chain table
%   bundled with scDock (receptor_reference.csv, human only), then the
%   AlphaFold model of the gene's reviewed UniProt entry, then the first
%   experimental structure UniProt lists for it.
%
%   The structure opens in MOLVIEWER, the Mol*-based viewer in the
%   Bioinformatics Toolbox (R2026b and later). When that is unavailable or
%   fails, the structure's page on RCSB or AlphaFold DB opens in the system
%   browser instead; both pages embed the same Mol* viewer. Local .pdb and
%   .pdbqt files fall back to viewPDBQT, which needs no toolbox.
%
%   INFO = gui.i_viewprotein(QUERY, Name=Value) sets options:
%     Viewer  - "auto" (default): molviewer, falling back as above
%               "molviewer": molviewer only; errors if it fails
%               "web": the RCSB / AlphaFold DB page in the system browser
%               "pdbqt": viewPDBQT from scDock; remote structures are
%               downloaded to tempdir first
%               "none": resolve the structure without opening anything
%     Species - "auto" (default), "human" or "mouse". "auto" treats an
%               all-uppercase symbol as human and anything else as mouse.
%     Source  - "auto" (default), "alphafold" or "pdb": which structure
%               to use for a gene. "alphafold" and "pdb" skip the curated
%               table.
%     ParentFig - figure to own a wait bar during the lookup and an error
%               dialog if it fails, for use from a GUI callback. With a
%               ParentFig, failures are reported in the dialog instead of
%               raised, and INFO.Viewer is "none". Default [] (no GUI).
%
%   INFO is a struct describing what was shown:
%     Query     - QUERY as a string
%     Kind      - "gene", "pdb" or "file"
%     Source    - "reference", "alphafold", "pdb" or "file"
%     Accession - UniProt accession ("" unless looked up)
%     Target    - PDB ID, structure URL, or file path handed to the viewer
%     Page      - web page for the structure ("" for a local file)
%     Viewer    - the viewer actually used ("molviewer", "web", "pdbqt"
%                 or "none")
%     Figure    - the viewer figure, or [] for "web" and "none"
%
%   Looking up a gene needs a network connection. An unresolvable gene
%   raises the error gui:i_viewprotein:NotFound.
%
%   Examples:
%     gui.i_viewprotein("EGFR")
%     gui.i_viewprotein("1IYT", Viewer="web")
%     gui.i_viewprotein("dock_output/receptor.pdbqt", Viewer="pdbqt")
%
%   See also molviewer, viewPDBQT, sc_dock_gene2pdb.

arguments
    query (1,1) string
    options.Viewer (1,1) string {mustBeMember(options.Viewer, ...
        ["auto", "molviewer", "web", "pdbqt", "none"])} = "auto"
    options.Species (1,1) string {mustBeMember(options.Species, ...
        ["auto", "human", "mouse"])} = "auto"
    options.Source (1,1) string {mustBeMember(options.Source, ...
        ["auto", "alphafold", "pdb"])} = "auto"
    options.ParentFig = []
end

parentfig = options.ParentFig;
if isempty(parentfig)
    info = i_run(query, options);
    return;
end
fw = gui.myWaitbar(parentfig);
try
    info = i_run(query, options);
    gui.myWaitbar(parentfig, fw);
catch ME
    gui.myWaitbar(parentfig, fw, true);
    gui.myErrordlg(parentfig, ME.message, "Protein Structure");
    info = i_emptyinfo(query);
end
end


function info = i_emptyinfo(query)
info = struct("Query", query, "Kind", "", "Source", "", ...
    "Accession", "", "Target", "", "Page", "", "Viewer", "none", ...
    "Figure", []);
end


function info = i_run(query, options)
% Resolve QUERY to a structure and show it; the body of i_viewprotein.
query = strtrim(query);
info = i_emptyinfo(query);

if isfile(query)
    info.Kind = "file";
    info.Source = "file";
    info.Target = query;
elseif i_ispdbid(query)
    info.Kind = "pdb";
    info.Source = "pdb";
    info.Target = upper(query);
    info.Page = i_rcsbpage(info.Target);
else
    info.Kind = "gene";
    info = i_resolvegene(info, options.Species, options.Source);
end

if options.Viewer == "none"
    return;
end
info = i_show(info, options.Viewer);
end


function info = i_resolvegene(info, species, source)
% Fill Source, Accession, Target and Page for a gene symbol.
gene = info.Query;
if species == "auto"
    if gene == upper(gene)
        species = "human";
    else
        species = "mouse";
    end
end

if source == "auto" && species == "human"
    pdbid = i_referencepdb(gene);
    if pdbid ~= ""
        info.Source = "reference";
        info.Target = pdbid;
        info.Page = i_rcsbpage(pdbid);
        return;
    end
end

[accession, pdbids] = i_uniprotlookup(gene, species);
if accession == ""
    error("gui:i_viewprotein:NotFound", ...
        "No reviewed %s UniProt entry was found for gene %s. " + ...
        "Check the gene symbol and the network connection.", species, gene);
end
info.Accession = accession;

if source ~= "pdb"
    modelUrl = i_alphafoldurl(accession);
    if modelUrl ~= ""
        info.Source = "alphafold";
        info.Target = modelUrl;
        info.Page = "https://alphafold.ebi.ac.uk/entry/" + accession;
        return;
    end
end
if ~isempty(pdbids)
    info.Source = "pdb";
    info.Target = pdbids(1);
    info.Page = i_rcsbpage(pdbids(1));
    return;
end
error("gui:i_viewprotein:NotFound", ...
    "UniProt entry %s (gene %s) has no AlphaFold model or PDB structure.", ...
    accession, gene);
end


function info = i_show(info, viewer)
% Open INFO.Target with VIEWER, falling back when VIEWER is "auto".
switch viewer
    case "molviewer"
        info.Figure = molviewer(i_molviewerinput(info.Target));
        info.Viewer = "molviewer";
    case "web"
        if info.Page == ""
            error("gui:i_viewprotein:NoPage", ...
                "A local file has no web page. Use Viewer=""auto"" or " + ...
                "Viewer=""pdbqt"" to view %s.", info.Target);
        end
        web(info.Page, "-browser");
        info.Viewer = "web";
    case "pdbqt"
        info.Figure = i_viewpdbqt(info.Target);
        info.Viewer = "pdbqt";
    case "auto"
        try
            info = i_show(info, "molviewer");
        catch ME
            warning("gui:i_viewprotein:MolviewerFailed", ...
                "molviewer could not open the structure (%s). " + ...
                "Using a fallback viewer instead.", ME.message);
            if info.Kind == "file"
                info = i_show(info, "pdbqt");
            else
                info = i_show(info, "web");
            end
        end
    otherwise
        % Unreachable: the arguments block restricts Viewer, and "none"
        % returns before i_show is called.
        error("gui:i_viewprotein:BadViewer", "Unknown viewer %s.", viewer);
end
end


function target = i_molviewerinput(target)
% molviewer reads PDB IDs, URLs and most structure formats, but not PDBQT.
% Convert a .pdbqt file to .pdb with OpenBabel, which scDock already needs.
if ~isfile(target)
    return;
end
[~, name, ext] = fileparts(target);
if lower(ext) ~= ".pdbqt"
    return;
end
pdbfile = fullfile(tempdir, name + ".pdb");
[status, msg] = system(sprintf('obabel "%s" -O "%s"', target, pdbfile));
if status ~= 0 || ~isfile(pdbfile)
    error("gui:i_viewprotein:ObabelFailed", ...
        "Converting %s to PDB with OpenBabel failed: %s", target, strtrim(msg));
end
target = string(pdbfile);
end


function fig = i_viewpdbqt(target)
% Show TARGET with scDock's viewPDBQT, downloading it first if remote.
if ~isfile(target)
    localfile = fullfile(tempdir, i_filename(target));
    websave(localfile, i_downloadurl(target));
    target = string(localfile);
end
if isempty(which("viewPDBQT"))
    run.ml_scDock();
end
[~, ~, ext] = fileparts(target);
if ~ismember(lower(ext), [".pdb", ".pdbqt"])
    error("gui:i_viewprotein:UnsupportedFile", ...
        "viewPDBQT reads only .pdb and .pdbqt files, not %s.", target);
end
viewPDBQT(char(target));
fig = gcf;
end


function url = i_downloadurl(target)
% A PDB ID becomes its RCSB download URL; a URL is used as it is.
if i_ispdbid(target)
    url = "https://files.rcsb.org/download/" + target + ".pdb";
else
    url = target;
end
end


function name = i_filename(target)
if i_ispdbid(target)
    name = target + ".pdb";
else
    [~, base, ext] = fileparts(target);
    name = base + ext;
end
end


function tf = i_ispdbid(s)
% A PDB ID is a digit followed by three alphanumerics, e.g. 1IYT.
tf = ~isempty(regexp(s, "^[1-9][A-Za-z0-9]{3}$", "once"));
end


function page = i_rcsbpage(pdbid)
page = "https://www.rcsb.org/3d-view/" + upper(pdbid);
end


function pdbid = i_referencepdb(gene)
% First PDB ID listed for GENE in scDock's receptor_reference.csv, or "".
pdbid = "";
reffile = fullfile(fileparts(fileparts(mfilename("fullpath"))), ...
    "external", "ml_scDock", "receptor_reference.csv");
if ~isfile(reffile)
    return;
end
T = readtable(reffile, TextType="string");
hit = T.PDB_model(T.protein_name == gene);
if isempty(hit) || ismissing(hit(1)) || hit(1) == ""
    return;
end
% Entries look like "1B6E.A; 2L35.A": keep the first, drop the chain.
pdbid = extractBefore(strtrim(extractBefore(hit(1) + ";", ";")) + ".", ".");
end


function [accession, pdbids] = i_uniprotlookup(gene, species)
% Reviewed UniProt accession for GENE and the PDB IDs it cross-references.
accession = "";
pdbids = strings(0, 1);
taxon = "9606";
if species == "mouse"
    taxon = "10090";
end
url = "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:" + ...
    gene + "+AND+organism_id:" + taxon + "+AND+reviewed:true" + ...
    "&fields=accession,xref_pdb&format=tsv&size=1";
try
    txt = string(webread(url, weboptions(Timeout=20)));
catch ME
    warning("gui:i_viewprotein:UniProtFailed", ...
        "UniProt lookup for %s failed: %s", gene, ME.message);
    return;
end
rows = splitlines(strtrim(txt));
if numel(rows) < 2
    return;
end
cols = split(rows(2), char(9));
accession = strtrim(cols(1));
if numel(cols) >= 2 && strtrim(cols(2)) ~= ""
    pdbids = strtrim(split(strip(cols(2), "right", ";"), ";"));
end
end


function url = i_alphafoldurl(accession)
% URL of the AlphaFold DB model file for ACCESSION, or "" if none. The API
% is asked rather than building the URL, because the model version in the
% file name changes between AlphaFold DB releases.
url = "";
try
    r = webread("https://alphafold.ebi.ac.uk/api/prediction/" + accession, ...
        weboptions(Timeout=20));
catch
    % No model (the API answers 404) or no network: the caller falls back
    % to an experimental structure.
    return;
end
if iscell(r)
    r = r{1};
end
if ~isempty(r) && isfield(r, "pdbUrl")
    url = string(r(1).pdbUrl);
end
end
