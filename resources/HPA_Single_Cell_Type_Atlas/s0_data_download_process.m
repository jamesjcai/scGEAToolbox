urltxt='https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip';
fn=unzip(urltxt);
T=readtable(fn{1},'FileType','text');
delete(fn{1});
