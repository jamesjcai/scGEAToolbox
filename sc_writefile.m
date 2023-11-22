function sc_writefile(filename, X, g, delim)
if nargin < 4
    delim = '\t';
end
t = table();
t.genes = string(g(:));
t = [t, array2table(X)];
%[file,path] = uiputfile('data_1.txt');
%if file~=0
writetable(t, filename, 'Delimiter', delim, 'filetype', 'text');
%end