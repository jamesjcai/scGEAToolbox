function [setmatrx, setnames, setgenes] = e_getenrichrsets(libraries, opts)
%E_GETENRICHRSETS  Load Enrichr gene set libraries as a membership matrix.
%
%   [setmatrx, setnames, setgenes] = E_GETENRICHRSETS() downloads the five
%   libraries the fGSEA workflow has always used and returns them in the
%   same form as PKG.E_GETGENESETS, ready for SC_GSETTEST.
%
%   [...] = E_GETENRICHRSETS(libraries) fetches the named libraries instead.
%   PKG.I_GET_ENRICHR_LIBRARIES lists what is available.
%
%   Downloads are cached on disk, because these files are a few megabytes
%   each and are refetched on every call otherwise. Set Refresh=true to
%   ignore the cache.
%
%   USAGE:
%     [M, names, genes] = pkg.e_getenrichrsets();
%     T = sc_gsettest(stats, genelist, M, names, genes, Method="gsea");
%
%   INPUT:
%     libraries - string array of Enrichr library names. Default is
%                 ["KEGG_2019_Human", "BioPlanet_2019",
%                  "GO_Biological_Process_2018",
%                  "GO_Molecular_Function_2018", "Reactome_2016"], the set
%                 the R fgsea script used, kept so results stay comparable.
%
%   NAME-VALUE:
%     CacheDir - where to keep the downloaded libraries (default a folder
%                under TEMPDIR).
%     Refresh  - re-download even when a cached copy exists (default false).
%     Prefix   - prepend the library name to each set name, so sets from
%                different libraries stay distinguishable (default true).
%     Verbose  - print progress (default false).
%
%   OUTPUTS:
%     setmatrx - nSets-by-nGenes sparse logical membership matrix.
%     setnames - nSets-by-1 set names.
%     setgenes - nGenes-by-1 gene symbols labelling the columns.
%
%   See also PKG.E_GETGENESETS, PKG.I_GET_ENRICHR_LIBRARIES, SC_GSETTEST.

arguments
    libraries string = ["KEGG_2019_Human", "BioPlanet_2019", ...
        "GO_Biological_Process_2018", "GO_Molecular_Function_2018", ...
        "Reactome_2016"]
    opts.CacheDir string = ""
    opts.Refresh (1,1) logical = false
    opts.Prefix (1,1) logical = true
    opts.Verbose (1,1) logical = false
end

if isempty(libraries)
    error("e_getenrichrsets:noLibraries", ...
        "Name at least one Enrichr library.");
end

cacheDir = opts.CacheDir;
if strlength(cacheDir) == 0
    cacheDir = fullfile(tempdir, "scgeatoolbox_enrichr");
end
if ~isfolder(cacheDir)
    mkdir(cacheDir);
end

allSets = {};
allNames = strings(0, 1);
for k = 1:numel(libraries)
    lines = i_fetchlibrary(libraries(k), cacheDir, opts.Refresh, opts.Verbose);
    [names, sets] = i_parsegmt(lines);
    if isempty(names)
        warning("e_getenrichrsets:emptyLibrary", ...
            "Library %s returned no gene sets and was skipped.", ...
            libraries(k));
        continue;
    end
    if opts.Prefix
        names = libraries(k) + ": " + names;
    end
    allNames = [allNames; names(:)]; %#ok<AGROW>
    allSets = [allSets; sets(:)]; %#ok<AGROW>
end

if isempty(allNames)
    error("e_getenrichrsets:noSets", ...
        "None of the requested libraries returned any gene sets.");
end

% One pass to build the universe, then one sparse matrix built from triplets
% rather than grown set by set.
setgenes = unique(vertcat(allSets{:}));
setgenes(strlength(setgenes) == 0) = [];
lookup = containers.Map(cellstr(setgenes), num2cell(1:numel(setgenes)));

rowIdx = cell(numel(allSets), 1);
colIdx = cell(numel(allSets), 1);
for k = 1:numel(allSets)
    g = allSets{k};
    g(strlength(g) == 0) = [];
    if isempty(g)
        rowIdx{k} = [];
        colIdx{k} = [];
        continue;
    end
    cols = cell2mat(values(lookup, cellstr(g)));
    colIdx{k} = cols(:);
    rowIdx{k} = repmat(k, numel(cols), 1);
end
rowIdx = vertcat(rowIdx{:});
colIdx = vertcat(colIdx{:});
setmatrx = sparse(rowIdx, colIdx, true, numel(allSets), numel(setgenes));
setnames = allNames;

if opts.Verbose
    fprintf("e_getenrichrsets: %d sets over %d genes from %d libraries\n", ...
        numel(setnames), numel(setgenes), numel(libraries));
end

end


function lines = i_fetchlibrary(library, cacheDir, refresh, verbose)
% Enrichr serves a library as tab-separated text, one gene set per line.
cacheFile = fullfile(cacheDir, matlab.lang.makeValidName(library) + ".gmt");
if ~refresh && isfile(cacheFile)
    if verbose
        fprintf("  %s (cached)\n", library);
    end
    lines = i_readlines(cacheFile);
    if ~isempty(lines)
        return;
    end
    % An empty cache file is a failed download from an earlier run; fall
    % through and fetch it again rather than returning nothing.
end

url = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text" + ...
    "&libraryName=" + library;
if verbose
    fprintf("  %s (downloading)\n", library);
end
options = weboptions(Timeout=120, ContentType="text");
text = webread(url, options);

% Write the cache only once the download has succeeded, so an interrupted
% run cannot leave a truncated file to be trusted next time.
tmpFile = cacheFile + ".part";
fid = fopen(tmpFile, "w", "n", "UTF-8");
if fid < 0
    error("e_getenrichrsets:cacheWrite", ...
        "Cannot write the cache file %s.", tmpFile);
end
closeFile = onCleanup(@() fclose(fid));
fwrite(fid, unicode2native(string(text), "UTF-8"));
clear closeFile;
movefile(tmpFile, cacheFile);

lines = i_readlines(cacheFile);
end


function lines = i_readlines(file)
lines = readlines(file, EmptyLineRule="skip");
end


function [names, sets] = i_parsegmt(lines)
% GMT: set name, then a description field Enrichr leaves empty, then the
% genes. Some libraries append a weight to each gene as "GENE,1.0".
names = strings(numel(lines), 1);
sets = cell(numel(lines), 1);
keep = false(numel(lines), 1);
for k = 1:numel(lines)
    parts = split(strip(lines(k), "right", newline), sprintf("\t"));
    if numel(parts) < 3
        continue;
    end
    genes = parts(3:end);
    genes = extractBefore(genes + ",", ",");
    genes = upper(strip(genes));
    genes(strlength(genes) == 0) = [];
    if isempty(genes)
        continue;
    end
    names(k) = strip(parts(1));
    sets{k} = unique(genes);
    keep(k) = true;
end
names = names(keep);
sets = sets(keep);
end
