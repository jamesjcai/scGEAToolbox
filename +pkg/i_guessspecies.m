function [a] = i_guessspecies(g)
%I_GUESSSPECIES Guess the species of a gene list from its symbol casing.
%   a = pkg.i_guessspecies(g) returns 'human' or 'mouse'. HGNC symbols are
%   upper case (ACTB, CD8A); MGI symbols are title case (Actb, Cd8a), so
%   the share of all-upper-case symbols separates the two cleanly. A few
%   upper-case symbols in a mouse list (MT-CO1 and friends) are tolerated
%   by the 0.9 threshold.
%
%   Accepts a string array, a cell array of character vectors, or a
%   character vector. The cellstr case used to error -- "==" on a cell
%   array is a numeric comparison -- which mattered once this function
%   acquired its first caller in gui.callback_SubtypeAnnotation.

g = string(g);
g = g(strlength(g) > 0);
if isempty(g)
    a = 'human';
    return
end

a = 'human';
if mean(upper(g) == g) < 0.9
    a = 'mouse';
end
end
