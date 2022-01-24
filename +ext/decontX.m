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
	if exist('./output.h5','file'), delete('./output.h5'); end
end
save('input.mat','X','-v7.3');
pkg.RunRcode('script.R');
if exist('./output.h5','file')
    X=h5read('output.h5','/X');
    contamination=h5read('output.h5','/contamination');
    % load('output.mat','X','contamination')
end
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.h5','file'), delete('./output.h5'); end
end
cd(oldpth);
end
