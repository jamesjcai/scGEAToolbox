function [T]=r_sctenifoldnet(sce1,sce2)

if ~isa(sce1,'SingleCellExperiment') || ~isa(sce2,'SingleCellExperiment')
    error('SCE1 and SCE2 should be SingleCellExperiment objects.');
end
T=[];
assert(isequal(sce1.g,sce2.g))
isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_scTenifoldNet');
if ~isok, error(msg); end

if ~isdebug
	if exist('./input1.rds','file'), delete('./input1.rds'); end
	if exist('./input2.rds','file'), delete('./input2.rds'); end
end

sc_sce2rds(sce1,fullfile(pwd(), 'input1.rds'));
sc_sce2rds(sce2,fullfile(pwd(), 'input2.rds'));

if ~isdebug
	if exist('./output.txt','file'), delete('./output.txt'); end
end
if isdebug, return; end
pkg.RunRcode('script.R');
if exist('./output.txt','file')
    T=readtable('output.txt');
end
if ~isdebug
	if exist('./input1.rds','file'), delete('./input1.rds'); end
	if exist('./input2.rds','file'), delete('./input2.rds'); end
	if exist('./output.txt','file'), delete('./output.txt'); end
end
cd(oldpth);
end
