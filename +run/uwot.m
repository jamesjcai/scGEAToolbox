function [s]=uwot(X)
%UWOT - R pacakge for UMAP
oldpth=pwd();
[isok,msg]=commoncheck_R('R_uwot');
if ~isok, error(msg); end
if issparse(X), X=full(X); end

if exist('output.csv','file'), delete('output.csv'); end
csvwrite('input.csv',X');

RunRcode('script.R');
if exist('output.csv','file')
    s=csvread('output.csv',1,1);
else
    s=[];
end
if exist('input.csv','file'), delete('input.csv'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end
