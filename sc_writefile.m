function sc_writefile(filename, X, g, delim)

if nargin < 4
    delim = '\t';
end
t = table();
t.genes = string(g(:));
t = [t, array2table(X)];
writetable(t, filename, 'Delimiter', delim, 'filetype', 'text');
