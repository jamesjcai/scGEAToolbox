function [X]=decontX(X)
%Run decontX 
%
% see also: soupX
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

if isa(X,'SingleCellExperiment')
    X=X.X;
else
    if nargin<2, error("[c]=run.decontX(X,genelist)"); end    
end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_decontX');
if ~isok, error(msg); end

if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end
writematrix(X,'input.txt');
RunRcode('script.R');

if exist('./output.csv','file')
    X=readmatrix('./output.csv');    
end
if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end
cd(oldpth);
end