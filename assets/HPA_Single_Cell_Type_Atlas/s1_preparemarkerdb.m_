Tx = T;
Tx(log(Tx.nTPM) < 2, :) = [];
a = string(unique(Tx.CellType));
gt = string(Tx.CellType);
genesymbollist = string(Tx.GeneName);

%%
fid = fopen('outfile.txt', 'w');
for k = 1:length(a)
    k;
    idx = find(gt == a(k));
    fprintf(fid, '%s\t', a(k));
    % fprintf(fid,'%s\t',organlist(k));
    fprintf(fid, '%s,', genesymbollist(idx));
    fprintf(fid, '\n');
end
fclose(fid);