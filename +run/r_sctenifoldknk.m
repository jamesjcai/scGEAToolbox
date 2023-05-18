function [T]=r_sctenifoldknk(sce,targetg)
T=[];
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
	if exist('./output.txt','file'), delete('./output.txt'); end
end
if exist('./input.h5','file'), delete('./input.h5'); end

h5create('input.h5', '/X', size(sce.X));
h5write('input.h5', '/X', sce.X);
h5create('input.h5', '/g', size(sce.g),'Datatype','string');
h5write('input.h5', '/g', sce.g);
h5create('input.h5', '/targetg', size(targetg),'Datatype','string');
h5write('input.h5', '/targetg', targetg);

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if exist('./output.txt','file')
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    T=readtable('output.txt');
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');
end
if ~isdebug
	if exist('./input.h5','file'), delete('./input.h5'); end
	if exist('./output.txt','file'), delete('./output.txt'); end
end
cd(oldpth);
end
