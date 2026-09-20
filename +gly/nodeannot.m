function T = nodeannot()
%NODEANNOT  Per-gene documented glycosylation site counts (UniProt).
%   T = GLY.NODEANNOT() loads the node-level glycosite table used by
%   GLY.LRWEIGHT to break ties within a ligand-receptor family: two
%   genes matched by the same Channel-B rule (e.g. FGFR1 and FGFR4, both
%   matched by the heparan-sulfate rule) previously received the identical
%   weight, because the rule matches by symbol PATTERN, not by gene. This
%   table gives each gene its own count, so the two can differ.
%
%   OUTPUT:
%     T - table, one row per gene in ASSETS/LIGAND_RECEPTOR/LIGAND_RECEPTOR.MAT
%         (ligand or receptor column, n=1347), with columns:
%           gene              - upper-case gene symbol
%           uniprot_accession - the reviewed (Swiss-Prot) human entry used
%           found             - true if a reviewed human UniProt entry was
%                               found for this symbol
%           n_nlinked         - count of CARBOHYD features whose note contains
%                               "N-linked"
%           n_o_galnac        - CARBOHYD, "O-linked ... GalNAc" (mucin-type)
%           n_o_fuc           - CARBOHYD, "O-linked ... Fuc" (as on Notch EGF
%                               repeats)
%           n_o_glcnac        - CARBOHYD, "O-linked ... GlcNAc"
%           n_o_other         - CARBOHYD, any other O-linked note (includes
%                               O-mannose and the O-xylose that primes GAG
%                               chain attachment, not split out further)
%           n_gpi             - 1 if a LIPID feature's note contains
%                               "GPI-anchor", else 0
%           n_total           - sum of the five site-count columns
%
%   WHAT A COUNT MEANS, AND WHAT IT DOES NOT. These are counts of sites
%   UniProt currently documents as sequence FEATURES (CARBOHYD/LIPID lines),
%   drawn from a mix of experimental evidence and sequence-motif prediction -
%   the raw TSV's evidence codes are not carried through to this table, so a
%   count does not distinguish the two. A gene with n_total=0 means no site is
%   currently ANNOTATED for it, not that the protein carries no glycan - most
%   membrane proteins are glycosylated somewhere, and UniProt's curation depth
%   varies by how well-studied the protein is. Read differences between two
%   well-studied family members (both heavily annotated) as more informative
%   than differences involving a sparsely-annotated one.
%
%   COVERAGE AND CURATION SHORTCUTS, STATED PLAINLY.
%     - 1346 of 1347 query genes matched a reviewed human UniProt entry via
%       gene:SYMBOL; CGB (chorionic gonadotropin beta) matched none, likely
%       because its UniProt primary symbol is a numbered paralog (CGB3/CGB5/
%       ...) with CGB itself only a synonym on some entries, not all.
%     - 58 of 1347 queries matched more than one reviewed entry (synonym
%       overlap between paralogs); the first returned entry was kept rather
%       than resolving to the best canonical isoform. Both of these are
%       simplifications a future pass could tighten; they were not judged
%       worth blocking this on.
%     - Data fetched from https://rest.uniprot.org/uniprotkb/search on the
%       date recorded in this file's SOURCE_NOTE variable. UniProt is updated
%       continuously; this is a snapshot, not a live query.
%
%   A gene not in the L-R database at all (not the receiver of a lookup this
%   function's caller would make) is simply absent from T; GLY.LRWEIGHT
%   treats a missing or unfound gene as zero documented sites, which in
%   CENTER mode is neutral only relative to whatever else that rule matched in
%   the same call - see its help for how the per-rule normalisation handles
%   this.
%
% see also: GLY.LRWEIGHT, GLY.LECTINMAP, ASSETS/GLYCOSYLATION

persistent cached

if isempty(cached)
    pw = fileparts(mfilename('fullpath'));
    matfile = fullfile(pw, '..', 'assets', 'Glycosylation', 'glycosite_annotations.mat');
    if ~isfile(matfile)
        error("GLY:NODEANNOT:AssetMissing", ...
            "%s not found. This asset ships with the toolbox; if it is " + ...
            "genuinely missing, re-run the fetch recorded in its source_note.", ...
            matfile);
    end
    S = load(matfile, 'T');
    cached = S.T;
end

T = cached;

end
