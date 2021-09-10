function [X,contamination]=decontX(X)
%Run decontX decontamination
%
% see also: run.SoupX
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_decontX');
if ~isok 
    error(msg); 
    X=[];
    contamination=[];
    return;
end

if isa(X,'SingleCellExperiment')
    X=X.X;
end

if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
end
save('input.mat','X','-v6');
pkg.RunRcode('script.R');
if exist('./output.mat','file')
    load('output.mat','X','contamination')
end
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
end
cd(oldpth);
end
