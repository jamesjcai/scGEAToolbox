[X,genelist]=sc_readfile('example_data\GSM3204304_P_P_Expr.csv');
[X,genelist]=sc_rmmtgenes(X,genelist);
[X,genelist]=sc_selectg(X,genelist,3,2);
[X]=sc_selectc(X,20000);
[X,genelist]=sc_selectg(X,genelist,3,2);
