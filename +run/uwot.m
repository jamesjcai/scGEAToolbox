function [s]=uwot(X)
%UWOT - R pacakge for UMAP
oldpth=pwd();
[isok,msg]=commoncheck_R('R_uwot');
if ~isok, error(msg); s=[]; return; end

if exist('output.csv','file'), delete('output.csv'); end
if issparse(X), X=full(X); end
writematrix(transpose(X),'input.csv');

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);


if exist('output.csv','file')
    s=readmatrix('output.csv');
else
    s=[];
end
if exist('input.csv','file'), delete('input.csv'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end

