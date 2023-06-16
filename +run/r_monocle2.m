function [t,s]=r_monocle2(X)
% Run Monocle pseudotime analysis
%[t_mono,s_mono]=run.monocle(X);

t=[]; s=[];
isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_monocle2');
if ~isok, error(msg); end

tmpfilelist={'input.h5','output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%save('input.mat','X','-v7.3');
pkg.e_writeh5(X,[],'input.h5');
Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);
if exist('./output.csv','file')
    dat=readmatrix('output.csv');
    t=dat(:,2);
    s=dat(:,3:end);
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end
