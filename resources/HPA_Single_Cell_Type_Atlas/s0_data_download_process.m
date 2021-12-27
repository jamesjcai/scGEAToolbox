% ref: https://www.proteinatlas.org/news/2021-07-28/a-single-cell-type-map-of-human-tissues
urltxt='https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip';
fn=unzip(urltxt);
T=readtable(fn{1},'FileType','text');
delete(fn{1});
