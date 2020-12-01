websave('TFome.xlsx','https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-020-0742-6/MediaObjects/41587_2020_742_MOESM3_ESM.xlsx');
t=readtable('TFome.xlsx','ReadVariableNames',false);
tfome_tfgenes=string(t.Var1);
% REF: https://www.nature.com/articles/s41587-020-0742-6
delete('TFome.xlsx');
clear t