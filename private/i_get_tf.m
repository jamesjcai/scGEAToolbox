websave('Homo_sapiens_TF.txt','http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF');
T=readtable('Homo_sapiens_TF.txt');
hs_tfgenes=string(T.Symbol);

websave('Mus_musculus_TF.txt','http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF');
T=readtable('Mus_musculus_TF.txt');
mm_tfgenes=string(T.Symbol);
