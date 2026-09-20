function key = i_typekey(s)
%I_TYPEKEY Comparison key for a cell type label.
%
%   key = pkg.i_typekey(s)
%
% Normalizes a cell type label so that the spellings that mean one type share
% one key: a bracketed citation dropped, lowercase, punctuation to spaces, the
% cluster index SCE.ASSIGNCELLTYPE appends dropped, the last word
% singularized, and a trailing "lymphocyte" read as "cell". So "T cells",
% "T-cell", "T cells_{3}" and "T lymphocytes" all key to "t cell".
%
% This is the key PKG.I_MATCHPRIMARYTYPE compares labels and primary types on,
% and the key PKG.I_SUBTYPEOVERLAP looks its table up by; it lives here so the
% two cannot drift apart.
%
% see also: pkg.i_normalizetypename, pkg.i_matchprimarytype,
%           pkg.i_subtypeoverlap

% A square-bracketed citation is provenance, not part of the name, and has to
% go before the rest runs rather than after. celltypes.xlsx spells five of its
% entries this way - "Neurons [PMID:27339989]", "Fibroblasts [PMID:32769974]".
% PKG.I_NORMALIZETYPENAME turns the brackets into spaces, which leaves "pmid"
% as the last word; the digits are then stripped as a cluster index, and the
% singularizer below works on "pmid" instead of "neurons". The key came out
% "neurons pmid", which matches no primary type by any rule, so a cluster
% annotated with one of those labels was invisible to subtype annotation.
%
% Square brackets only. The subtype formats SC_CSUBTYPEANNO writes are
% "T cells (Regulatory)" and "T cells_{Regulatory}", and a rule that ate a
% trailing parenthesis would key the first of them to a bare "t cell" - that
% is, would report an already-subtyped cell as one still waiting to be
% subdivided.
key = regexprep(string(s), '\s*\[[^\]]*\]\s*$', '');

key = pkg.i_normalizetypename(key);
key = regexprep(key, '\s+\d+$', '');            % 'T cells_{3}' -> 't cells'
key = regexprep(key, '(\w)s$', '$1');           % 't cells' -> 't cell'
key = regexprep(key, '\<lymphocyte$', 'cell');  % 'B lymphocytes' -> 'b cell'
key = strtrim(key);
end
