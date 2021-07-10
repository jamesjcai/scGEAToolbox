function [s]=uwot(X)
%UWOT - R pacakge for UMAP
oldpth=pwd();
[isok,msg]=commoncheck_R('R_uwot');
if ~isok, error(msg); end

if exist('output.csv','file'), delete('output.csv'); end
writematrix(transpose(X),'input.csv');

pkg.RunRcode('script.R');
if exist('output.csv','file')
    s=readmatrix('output.csv');
else
    s=[];
end
if exist('input.csv','file'), delete('input.csv'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end

