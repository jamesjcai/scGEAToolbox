function [X2]=r_SuperCell(X,gammavalue,kvalue)
% SuperCell  - gene expression recovery for single-cell RNA sequencing
% https://github.com/GfellerLab/SuperCell

if nargin<3, kvalue=5; end        % number of nearest neighbors to build kNN network
if nargin<2, gammavalue=5; end    % number of cells merge (grianing level)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SuperCell');
if ~isok, error(msg); X2=[]; return; end

if exist('./output.mat.tmp','file'), delete('./output.mat.tmp'); end

if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
end
lastwarn('');
save('input.mat','X','gammavalue','kvalue','-v6');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)    
    disp(warnId)
    error(warnId);
end

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if exist('output.mat','file')
    load('output.mat','X2');
end
if ~isdebug
	if exist('input.mat','file'), delete('input.mat'); end
	if exist('output.mat','file'), delete('output.mat'); end
end
cd(oldpth);
end