function s = i_normalizetypename(s)
% I_NORMALIZETYPENAME - Normalize a label for fuzzy cell type name comparison.
%
%   s = pkg.i_normalizetypename(s)
%
% Lowercases the input and replaces every run of non-alphanumeric characters
% with a single space, so that "CD8+ T-cells", "CD8 T cells" and "cd8_t_cells"
% all collapse to the same key.
%
% see also: pkg.i_celltypevocab, pkg.i_guesscelltypecol

s = lower(string(s));
s = regexprep(s, '[^a-z0-9]+', ' ');
s = strtrim(s);
end
