function [output]=i_ensembl2symbol(input)

if nargin<1
    input=["ENSG00000121410","ENSG00000121411"];
end

a=tempname;
b=websave(a,'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt');
T=readtable(b,'FileType','text');
% string([T.ensembl_gene_id, T.symbol, T.name])
keySet=T.ensembl_gene_id;
valueSet=T.symbol;
output=input;

if isstring(input)
    M = containers.Map(string(keySet),string(valueSet));
    %input=string({'ENSG00000121410', 'ENSG00000121411'});
    ix=isKey(M, cellstr(input));
    output(ix)=values(M,cellstr(input(ix)));
elseif iscell(input)
    disp('cell')
    M = containers.Map(keySet,valueSet);
    %input=string({'ENSG00000121410', 'ENSG00000121411'});
    ix=isKey(M, input);
    output(ix)=values(M,input(ix));    
end

