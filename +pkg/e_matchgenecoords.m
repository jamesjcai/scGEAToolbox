function [tf, loc, T] = e_matchgenecoords(geneNames, genome)
%E_MATCHGENECOORDS Match gene identifiers to genomic coordinate rows.
%
%   [tf, loc, T] = PKG.E_MATCHGENECOORDS(geneNames, genome) matches
%   GENENAMES, which may hold gene symbols or Ensembl gene IDs, against the
%   coordinate table T for GENOME. TF marks the entries that were placed
%   and LOC gives their row in T, exactly as ISMEMBER returns them.
%
%   Matching is case sensitive by design. Mouse and human symbols are
%   largely the same letters in different case -- Xkr4 against XKR4 -- so
%   folding case maps 54% of the mouse gene set onto human coordinates.
%   That is far too many to look like a failure and produces a confident,
%   entirely fictitious result. Case is the only thing separating the two
%   namespaces, so it is preserved.
%
%   For the same reason a symbol list whose casing contradicts the
%   requested genome is rejected outright rather than matched partially:
%   98.8% of human symbols are all-uppercase against 0.6% of mouse ones, so
%   the convention identifies the species almost perfectly.
%
%   If case-sensitive matching places under a quarter of the genes and
%   case-insensitive matching would do better, the latter is used and a
%   warning is issued. That covers a dataset which has upper-cased its own
%   symbols, without opening the cross-species door.
%
%   See also SC_INFERCNV, PKG.E_GENECOORDINATES, RUN.ML_SCEVAN.

arguments
    geneNames (:, 1) string
    genome (1, :) char
end

T = pkg.e_genecoordinates(genome);

if mean(startsWith(geneNames, "ENS", 'IgnoreCase', true)) > 0.5
    % ENSG versus ENSMUSG already separates the species unambiguously.
    key = T.EnsemblID;
else
    key = T.Gene;
    i_checkspecies(geneNames, genome);
end

[tf, loc] = ismember(geneNames, key);

if mean(tf) < 0.25
    [tfLoose, locLoose] = ismember(upper(geneNames), upper(key));
    if sum(tfLoose) > sum(tf)
        warning('pkg:e_matchgenecoords:caseFolded', ...
            ['Only %d of %d genes matched %s exactly; case-insensitive ' ...
            'matching places %d. Using it, but check that GENENAMES holds ' ...
            'identifiers of the expected species.'], ...
            sum(tf), numel(geneNames), genome, sum(tfLoose));
        tf = tfLoose;
        loc = locLoose;
    end
end

if ~any(tf)
    error('pkg:e_matchgenecoords:noMatch', ...
        ['None of the %d genes carry coordinates in genome %s. Check that ' ...
        'GENENAMES holds gene symbols or Ensembl IDs of the right species.'], ...
        numel(geneNames), genome);
end
end


function i_checkspecies(geneNames, genome)
% Human symbols are upper case, mouse symbols are capitalised. Refuse the
% combination that would otherwise half-match its way to a wrong answer.
hasLetters = geneNames ~= upper(geneNames) | geneNames ~= lower(geneNames);
symbols = geneNames(hasLetters);
if numel(symbols) < 20
    return;  % too few to read a convention from
end

looksHuman = mean(symbols == upper(symbols)) > 0.5;
wantsHuman = ~strcmp(genome, 'mm10');
if looksHuman == wantsHuman
    return;
end

if wantsHuman
    error('pkg:e_matchgenecoords:speciesMismatch', ...
        ['GENENAMES looks like mouse symbols but GENOME is %s. Pass ' ...
        'Genome="mm10" for mouse data, or convert the symbols to human ' ...
        'orthologues first.'], genome);
else
    error('pkg:e_matchgenecoords:speciesMismatch', ...
        ['GENENAMES looks like human symbols but GENOME is mm10. Pass ' ...
        'Genome="hg38" or Genome="hg19" for human data.']);
end
end
