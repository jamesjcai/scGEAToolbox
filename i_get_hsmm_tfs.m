warning off
websave('Homo_sapiens_TF.txt','http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF');
t=readtable('Homo_sapiens_TF.txt');
hs_tfgenes=string(t.Symbol);
%delete('Homo_sapiens_TF.txt');

websave('Mus_musculus_TF.txt','http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF');
t=readtable('Mus_musculus_TF.txt');
mm_tfgenes=string(t.Symbol);
%delete('Mus_musculus_TF.txt');
warning on
clear t