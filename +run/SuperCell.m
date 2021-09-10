function [X2]=SuperCell(X)
% SuperCell  - gene expression recovery for single-cell RNA sequencing
% https://github.com/GfellerLab/SuperCell

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SuperCell');
if ~isok, error(msg); X2=[]; return; end

if ~isdebug
	if exist('input.mat','file'), delete('input.mat'); end
	if exist('output.mat','file'), delete('output.mat'); end
end
save('input.mat','X');
pkg.RunRcode('script.R');
if exist('output.mat','file')
    load('output.mat','X2');
end
if ~isdebug
	if exist('input.mat','file'), delete('input.mat'); end
	if exist('output.mat','file'), delete('output.mat'); end
end
cd(oldpth);
end