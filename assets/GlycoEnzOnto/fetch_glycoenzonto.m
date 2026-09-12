function fetch_glycoenzonto()
%FETCH_GLYCOENZONTO  Refresh GlycoEnzOnto.gmt from its upstream repository.
%   Downloads the GMT and the licence text from github.com/neel-lab/
%   GlycoEnzOnto and writes both next to this file, which is where
%   GLY.ENZONTO reads them. Requires internet access.
%
%   The shipped copy is a dated snapshot. The upstream repository has been
%   quiet since November 2022, so this is unlikely to change under you -
%   but if it does, GLY.ENZONTO derives its aggregate flags and class
%   labels FROM THE GENE SETS rather than from a stored table, so a refresh
%   propagates without any further step. Run GLY.ENZONTO and check
%   the counts afterwards regardless: 122 terms, 27 of them aggregates, 403
%   genes, 95 leaf pathways over 370 genes as of the shipped snapshot.
%
%   LICENCE. GlycoEnzOnto is CC BY 4.0. It is redistributed inside this
%   toolbox with attribution to Groth, Diehl, Gunawan and Neelamegham,
%   "GlycoEnzOnto: a GlycoEnzyme pathway and molecular function ontology",
%   Bioinformatics 38 (2022) 5413-5420. Do not replace LICENSE.md with
%   anything but the upstream text.
%
%   Not called by anything else in the toolbox.
%
% see also: GLY.ENZONTO, GLY.GENESETS

here = fileparts(mfilename('fullpath'));
base = "https://raw.githubusercontent.com/neel-lab/GlycoEnzOnto/main/";
files = ["GlycoEnzOnto.gmt", "LICENSE.md"];

for k = 1:numel(files)
    dest = fullfile(here, files(k));
    fprintf('fetching %s ... ', files(k));
    try
        websave(dest, base + files(k));
    catch ME
        fprintf('FAILED\n');
        error("FETCH_GLYCOENZONTO:Download", ...
            "Could not fetch %s: %s", files(k), ME.message);
    end
    d = dir(dest);
    fprintf('%.0f KB\n', d.bytes / 1024);
end

% Fail loudly on an empty or obviously wrong GMT rather than leaving the
% asset broken for GLY.ENZONTO to trip over later.
lines = readlines(fullfile(here, "GlycoEnzOnto.gmt"));
lines = lines(strlength(strtrim(lines)) > 0);
if numel(lines) < 50
    error("FETCH_GLYCOENZONTO:Suspect", ...
        "GlycoEnzOnto.gmt has only %d non-empty lines; expected ~122.", ...
        numel(lines));
end
fprintf('\n%d pathway terms written. Verifying the parse ...\n', numel(lines));

clear gly.enzonto   % drop the persistent cache
[~, leafNames, leafGenes, T] = gly.enzonto();
[~, ~, allGenes, TA] = gly.enzonto(Include = "all");
fprintf('  %d terms, %d aggregates, %d genes\n', height(TA), ...
    sum(TA.IsAggregate), numel(allGenes));
fprintf('  %d leaf pathways over %d genes\n', numel(leafNames), numel(leafGenes));
fprintf('  classes: %s\n', strjoin(unique(T.Class(strlength(T.Class) > 0))', ', '));

end
