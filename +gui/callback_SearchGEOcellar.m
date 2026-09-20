function [sce] = callback_SearchGEOcellar(src, ~)
% CALLBACK_SEARCHGEOCELLAR  Search the GEOcellar index and import the result.
%
%   sce = gui.callback_SearchGEOcellar(src)
%
%   Called from the scGEAToolbox main menu (File -> Search GEOcellar to
%   Import Data...). Asks for a research topic, queries the GEOcellar
%   ChromaDB index of GEO metadata, and imports whichever hits are chosen.
%
%   The index holds two kinds of record and most hits are the second kind:
%     * a study (GSE), whose samples are all imported and merged
%     * a sample (GSM), which is imported on its own, fetched from the
%       series named in its ParentStudy metadata
%   Either way C_BATCH_ID holds the GSM accession each cell came from.
%
%   More than one result can be picked. The picks are grouped by series, so
%   several samples of one study are fetched as one study, and the studies
%   are then merged on their common genes into a single SCE. Picking a
%   study together with some of its own samples still asks for the study
%   once.
%
%   Requires passkey validation (Help → Validate Passkey); returns []
%   with an explanatory dialog when the passkey is missing.
%
%   Returns [] when the user cancels at any step, or when the search or the
%   import fails; the caller leaves the current dataset alone in that case.
%
%   The data comes from the GEOcellar GCS bucket, which holds a curated
%   subset of GEO. A record found in the index may name samples that were
%   never deposited there; those are skipped and reported rather than
%   failing the import.
%
%   See also LLM.GEOCELLAR.CHROMA_QUERY, LLM.GEOCELLAR.IMPORT_STUDY,
%   GUI.CALLBACK_LAUNCHGEOCELLAR

sce = [];
[FigureHandle, ~] = gui.gui_getfigsce(src);

if ~pkg.i_license
    gui.myErrordlg(FigureHandle, ...
        "This function requires passkey validation. You can " + ...
        "validate your passkey by selecting Help → Validate " + ...
        "Passkey from the menu.");
    return;
end

cfg = llm.geocellar.geocellar_config();
if isempty(char(cfg.ChromaAPIKey))
    gui.myErrordlg(FigureHandle, ['No ChromaDB API key found. Searching ' ...
        'GEOcellar needs a CHROMA_API_KEY entry in your llm_api_key.env ' ...
        'file, or a CHROMA_API_KEY environment variable set before ' ...
        'MATLAB started. Add the key, then try again.'], ...
        'gui:callback_SearchGEOcellar:noApiKey');
    return;
end

preftagname = 'geocellarsearchtopic';
previoustopic = getpref('scgeatoolbox', preftagname, 'T cell exhaustion');

answer = gui.myInputdlg({'Research topic:'}, 'Search GEOcellar', ...
    {previoustopic}, FigureHandle);
if isempty(answer), return; end
topic = strtrim(answer{1});
if isempty(topic), return; end
setpref('scgeatoolbox', preftagname, topic);

fw = gui.myWaitbar(FigureHandle, [], false, ...
    sprintf('Searching GEOcellar for "%s"...', topic));
try
    hits = llm.geocellar.chroma_query(topic, cfg, 20);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(hits)
    gui.myHelpdlg(FigureHandle, ...
        sprintf('Nothing in the GEOcellar index matched "%s".', topic));
    return;
end

items = i_describehits(hits);
prompt = sprintf(['%d results match "%s". Picking a study (GSE) imports ' ...
    'all of its samples merged into one dataset; picking a sample (GSM) ' ...
    'imports that sample alone. Ctrl-click or shift-click to pick several ' ...
    '- they are merged on their common genes. The GSM accession of each ' ...
    'cell is kept in SCE.C_BATCH_ID.'], numel(hits), topic);
ctxitems = struct( ...
    "Label", {@(k) i_geolabel(hits, k)}, ...
    "Callback", {@(k) i_opengeopage(hits, k)});
[indx, tf] = gui.myListdlg(FigureHandle, items, 'GEOcellar Results', ...
    [], true, true, [640, 420], prompt, true, ctxitems);
if tf ~= 1, return; end

[targets, bad] = i_resolvetargets(hits(indx));
if isempty(targets)
    gui.myErrordlg(FigureHandle, sprintf(['%s cannot be imported. A ' ...
        'sample record needs a ParentStudy in its metadata to say which ' ...
        'series to fetch it from, and none of these has one. Pick another ' ...
        'result, or import it by accession with File > Import Data.'], ...
        strjoin(bad, ', ')), ...
        'gui:callback_SearchGEOcellar:noParentStudy');
    return;
end

[sces, infos, failed, reasons] = i_importall(FigureHandle, targets, ...
    string(cfg.DataDir));

if isempty(sces)
    gui.myErrordlg(FigureHandle, ...
        char(strjoin(failed+": "+reasons, [newline newline])), ...
        'gui:callback_SearchGEOcellar:noneImported');
    return;
end

try
    sce = i_mergestudies(sces);
catch ME
    sce = [];
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

i_reportgaps(FigureHandle, infos, failed, reasons);
end


function label = i_geolabel(hits, k)
% I_GEOLABEL Say which records the context menu entry will open.
%   One accession is named outright, several are counted - a label listing
%   six accessions is wider than the dialog. Either way the label states the
%   number of tabs about to appear, since that is the part worth knowing
%   before clicking. "" hides the entry, which is what a selection holding
%   no openable record gets.

acc = i_accessions(hits, k);
if isempty(acc)
    label = "";
elseif isscalar(acc)
    label = "View " + acc + " on the GEO website";
else
    label = sprintf("View %d selected records on the GEO website", numel(acc));
end
end


function i_opengeopage(hits, k)
% I_OPENGEOPAGE Open each selected record's GEO page in the system browser.
%   The system browser, not MATLAB's: acc.cgi is CAPTCHA-gated, and the
%   built-in browser trips it where an ordinary browser with the user's own
%   cookies does not.

acc = i_accessions(hits, k);
for j = 1:numel(acc)
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + acc(j);
    try
        web(url, "-browser");
    catch
        % No system browser could be launched; the URL is still worth
        % handing over, and the rest of the selection still deserves a try
        try
            web(url);
        catch ME
            warning("gui:callback_SearchGEOcellar:noBrowser", ...
                "Could not open %s: %s", url, ME.message);
        end
    end
end
end


function acc = i_accessions(hits, k)
% I_ACCESSIONS The accessions of the rows an action was invoked on.
%   A record with no accession is dropped rather than opened as a truncated
%   URL, so the count I_GEOLABEL shows is the number of pages that will
%   actually open.

acc = strings(0, 1);
if isempty(k), return; end
k = k(k >= 1 & k <= numel(hits));
for j = 1:numel(k)
    one = i_text(hits(k(j)), "id");
    if strlength(one) > 0
        acc(end+1) = one; %#ok<AGROW>
    end
end
end


function [sces, infos, failed, reasons] = i_importall(FigureHandle, targets, dataDir)
% I_IMPORTALL Import each selected study, one waitbar across the lot.
%   A study that cannot be fetched is recorded and skipped rather than
%   losing the ones that can, the same bargain IMPORT_STUDY makes for an
%   individual sample. FAILED and REASONS pair up for the report.

sces = cell(numel(targets), 1);
infos = cell(numel(targets), 1);
failed = strings(0, 1);
reasons = strings(0, 1);

fw = gui.myWaitbar(FigureHandle, [], false, ...
    sprintf('Downloading %s from GEOcellar...', targets(1).id));
for k = 1:numel(targets)
    if k > 1
        fw = gui.myWaitbar(FigureHandle, fw, false, '', ...
            sprintf('Downloading %s from GEOcellar (%d of %d)...', ...
            targets(k).id, k, numel(targets)));
    end
    try
        [sces{k}, infos{k}] = llm.geocellar.import_study(targets(k), dataDir);
    catch ME
        failed(end+1) = string(targets(k).id); %#ok<AGROW>
        reasons(end+1) = string(ME.message); %#ok<AGROW>
    end
end
gui.myWaitbar(FigureHandle, fw);

keep = ~cellfun(@isempty, sces);
sces = sces(keep);
infos = infos(keep);
end


function sce = i_mergestudies(sces)
% I_MERGESTUDIES Merge one SCE per study, keeping the GSM batch IDs.
%   SC_MERGESCES rewrites C_BATCH_ID two ways: it renumbers 1..N when the
%   merged set holds one unique value, and it appends a _k suffix when two
%   inputs share a value. Neither should fire here, since the targets are
%   distinct series and a GSM belongs to one of them - so the accessions
%   are collected first and put back after, which is a no-op in the
%   expected case and keeps them readable if SC_MERGESCES ever changes its
%   mind. Cells are concatenated in input order, so they line up.

if isscalar(sces)
    sce = sces{1};
    return;
end

ids = cellfun(@(s) string(s.c_batch_id(:)), sces, "UniformOutput", false);
batchid = vertcat(ids{:});

sce = sc_mergesces(sces, "intersect", true);

if numel(sce.c_batch_id) == numel(batchid)
    sce.c_batch_id = batchid;
else
    warning("gui:callback_SearchGEOcellar:batchIdLost", ...
        ['Merged cell count does not match the imported samples; ' ...
        'leaving SCE.C_BATCH_ID as SC_MERGESCES wrote it.']);
end
end


function i_reportgaps(FigureHandle, infos, failed, reasons)
% I_REPORTGAPS One dialog for everything the selection did not get.
%   Samples absent from the bucket and whole studies that failed are the
%   same kind of news to the user, so they share a dialog rather than
%   queueing two.

lines = strings(0, 1);

for k = 1:numel(infos)
    info = infos{k};
    if isempty(info.skipped), continue; end
    nTotal = numel(info.loaded) + numel(info.skipped);
    lines(end+1) = sprintf(['%s: imported %d of %d samples. Not ' ...
        'available in the GEOcellar bucket: %s'], info.gse, ...
        numel(info.loaded), nTotal, ...
        strjoin(string({info.skipped.gsm}), ', ')); %#ok<AGROW>
end

for k = 1:numel(failed)
    lines(end+1) = sprintf('%s: not imported. %s', failed(k), reasons(k)); %#ok<AGROW>
end

if isempty(lines), return; end

gui.myHelpdlg(FigureHandle, char(strjoin(lines, [newline newline])), ...
    'Some Samples Skipped');
end


function [targets, bad] = i_resolvetargets(hits)
% I_RESOLVETARGETS Turn the selected hits into IMPORT_STUDY arguments.
%   Hits are grouped by series, so picking three samples of one study is
%   one pass over that study rather than three, and merging them is left to
%   IMPORT_STUDY. A study hit contributes its whole sample list, a sample
%   hit contributes only itself; IMPORT_STUDY parses the accessions out of
%   the joined text and drops duplicates, so a study picked alongside its
%   own samples still resolves to one request for the study. BAD lists the
%   hits that are neither - a sample with no ParentStudy names nothing to
%   fetch it from.

targets = struct("id", {}, "samples", {});
bad = strings(0, 1);
gseList = strings(0, 1);
parts = strings(0, 1);

for k = 1:numel(hits)
    acc = i_text(hits(k), "id");
    parent = i_text(hits(k), "parent_study");

    if startsWith(acc, "GSE")
        gse = acc;
        samples = i_text(hits(k), "samples");
    elseif startsWith(acc, "GSM") && startsWith(parent, "GSE")
        gse = parent;
        samples = acc;
    else
        bad(end+1) = acc; %#ok<AGROW>
        continue;
    end

    where = find(gseList == gse, 1);
    if isempty(where)
        gseList(end+1) = gse; %#ok<AGROW>
        parts(end+1) = samples; %#ok<AGROW>
    else
        parts(where) = parts(where) + "," + samples;
    end
end

for k = 1:numel(gseList)
    targets(k) = struct("id", char(gseList(k)), "samples", char(parts(k)));
end
end


function items = i_describehits(hits)
% I_DESCRIBEHITS One list line per hit: accession, what it is, organism, title.
items = strings(numel(hits), 1);
for k = 1:numel(hits)
    acc = i_text(hits(k), "id");
    ttl = i_text(hits(k), "title");
    org = i_text(hits(k), "organism");

    line = acc + " (" + i_kind(hits(k)) + ")";
    if strlength(org) > 0
        line = line + " [" + org + "]";
    end
    if strlength(ttl) > 0
        line = line + " " + ttl;
    end
    items(k) = line;
end
end


function kind = i_kind(hit)
% I_KIND How the list describes a hit: a study of N samples, or one sample.
acc = i_text(hit, "id");
if startsWith(acc, "GSE")
    n = i_countsamples(hit);
    if n > 0
        kind = sprintf("study, %d sample%s", n, string(repmat('s', 1, n ~= 1)));
    else
        kind = "study";
    end
    return;
end

parent = i_text(hit, "parent_study");
if strlength(parent) > 0
    kind = "sample of " + parent;
else
    kind = "sample";
end
end


function txt = i_text(hit, fieldname)
% I_TEXT Read one metadata field, tolerating its absence or emptiness.
%   The collection has entries from more than one ingest script, so a hit
%   is not guaranteed to carry every field.
txt = "";
if isfield(hit, fieldname) && ~isempty(hit.(fieldname))
    txt = strtrim(string(hit.(fieldname)));
    txt = strjoin(splitlines(txt), " ");
end
end


function n = i_countsamples(hit)
% I_COUNTSAMPLES How many GSM accessions the samples field lists.
%   Counted the way LLM.GEOCELLAR.IMPORT_STUDY parses them, so the number
%   shown is the number that will be attempted.
n = numel(unique(string(regexp(i_text(hit, "samples"), ...
    "GSM\d+", "match")), "stable"));
end
