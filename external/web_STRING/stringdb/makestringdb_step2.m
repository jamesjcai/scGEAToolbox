% T=readtable('stringdb.txt');
if ~exist('T', 'var')
    if isfile('stringdb_mouse.txt')
        T = readtable('stringdb_mouse.txt');
    else
        error('Input table T is missing and stringdb_mouse.txt was not found.');
    end
end

G = graph(T.protein1, T.protein2, T.combined_score);
% g=string(table2array(G.Nodes));
% G=G.rmnode(2);
save stringdb_mouse G