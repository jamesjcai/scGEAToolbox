function [X,contamination]=decontX(X)
%Run decontX 
%
% see also: soupX
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

if isa(X,'SingleCellExperiment')
    X=X.X;
end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_decontX');
if ~isok, error(msg); end
if exist('./input.mat','file'), delete('./input.mat'); end
if exist('./output.mat','file'), delete('./output.mat'); end
save('input.mat','X','-v6');
pkg.RunRcode('script.R');
if exist('./output.mat','file')
    load('output.mat','X','contamination')
end
if exist('./input.mat','file'), delete('./input.mat'); end
if exist('./output.mat','file'), delete('./output.mat'); end
cd(oldpth);
end
