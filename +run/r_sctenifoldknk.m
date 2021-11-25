function [res]=r_sctenifoldknk(sce,targetg)
res=[];
if nargin<2
    error('Need two input variables.');
end
assert(ismember(targetg,sce.g));

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_scTenifoldKnk');
if ~isok, error(msg); end
if ~isa(sce,'SingleCellExperiment')
    error('SCE should be a SingleCellExperiment object.');
end

if ~isdebug
	if exist('./input.h5','file'), delete('./input.h5'); end
	if exist('./output.h5','file'), delete('./output.h5'); end
end
if exist('./input.h5','file'), delete('./input.h5'); end

h5create('input.h5', '/X', size(sce.X));
h5write('input.h5', '/X', sce.X);
h5create('input.h5', '/g', size(sce.g),'Datatype','string');
h5write('input.h5', '/g', sce.g);
h5create('input.h5', '/targetg', size(targetg),'Datatype','string');
h5write('input.h5', '/targetg', targetg);
pkg.RunRcode('script.R');
if exist('./output.h5','file')
% --------
end
if ~isdebug
	if exist('./input.h5','file'), delete('./input.h5'); end
	if exist('./output.h5','file'), delete('./output.h5'); end
end
cd(oldpth);
end