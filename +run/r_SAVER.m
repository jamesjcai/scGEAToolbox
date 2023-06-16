function [X2]=r_SAVER(X)
% SAVER - gene expression recovery for single-cell RNA sequencing
% https://mohuangx.github.io/SAVER/articles/saver-tutorial.html

oldpth=pwd();
[isok,msg]=commoncheck_R('R_SAVER');
if ~isok, error(msg); X2=[]; return; end
if exist('input.mat','file'), delete('input.mat'); end
if exist('output.mat','file'), delete('output.mat'); end
save('input.mat','X');
Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);
if exist('output.mat','file')
    load('output.mat','X2');
end
if exist('input.mat','file'), delete('input.mat'); end
if exist('output.mat','file'), delete('output.mat'); end
cd(oldpth);
end