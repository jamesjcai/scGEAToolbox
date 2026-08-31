function X = i_extractgenerows(names, gU, Xside)
%I_EXTRACTGENEROWS  Rows of Xside matching upper-case gene symbols NAMES,
%in NAMES order. Xside is genes-by-cells for the population NAMES was
%selected from (Xs for sender variables, Xr for receiver variables).
[~, idx] = ismember(upper(names(:)), gU);
X = Xside(idx, :);
end
