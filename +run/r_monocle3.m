function [t]=r_monocle3(X,idx)
% Run Monocle pseudotime analysis
%[t]=run.r_monocle3(X,idx);
if nargin<2, idx=[1 2]; end
t=[];
isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_monocle3');
if ~isok, error(msg); end

tmpfilelist={'input.h5','output.h5','input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

save('input.mat','X','idx','-v7.3');
%pkg.e_writeh5(X,[],'input.h5');
Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);
if exist('./output.h5','file')
    t=h5read('output.h5','/t');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
