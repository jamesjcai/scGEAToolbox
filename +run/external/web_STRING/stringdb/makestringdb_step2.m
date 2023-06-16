% T=readtable('stringdb.txt');
G=graph(T.protein1,T.protein2,T.combined_score);
% g=string(table2array(G.Nodes));
% G=G.rmnode(2);
save stringdb_mouse G