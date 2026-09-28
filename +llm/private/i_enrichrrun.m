function [Tlist, err] = i_enrichrrun(genes, background, genesets, warnId)
%I_ENRICHRRUN Run Enrichr on one gene list; report an API failure, don't hide it.
%   [Tlist, err] = i_enrichrrun(genes, background, genesets, warnId) returns
%   one table per gene-set library, and ERR, the API's error message ('' on
%   success). An outage used to return the same empty tables as "no enriched
%   terms", with only a console warning, so the report read as a negative
%   result. Shared by llm.run_enrichr, llm.run_enrichr_sctenifoldnet and
%   llm.run_enrichr_sctenifoldknk, which each carried a copy.
Tlist = cell(numel(genesets), 1);
err = '';
if isempty(genes)
    return;
end
try
    Tlist = run.ml_Enrichr(genes, background, genesets);
catch ME
    err = ME.message;
    warning(warnId, 'Enrichr API error: %s', ME.message);
end
end
