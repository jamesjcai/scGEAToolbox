function [g]=get_ribosomalgenes
%Get ribosomal genes

options = weboptions('Timeout',21);
fname=tempname;
websave(fname,'https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch',options);
% t=readtable('a.txt','PreserveVariableNames',true);
warning off 
t=readtable(fname,'Range','A:B');
warning on
g=string(t.ApprovedSymbol);
delete(fname);
