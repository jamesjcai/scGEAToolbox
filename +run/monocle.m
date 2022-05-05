function [t,s]=monocle(X)
% Run Monocle pseudotime analysis
%[t_mono,s_mono]=run.monocle(X);

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_monocle');
if ~isok, error(msg); t=[]; s=[]; return; end

tmpfilelist={'input.mat','output.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

save('input.mat','X','-v7.3');
pkg.RunRcode('script.R');
if exist('./output.csv','file')
    dat=readmatrix('output.csv');
    t=dat(:,2);
    s=dat(:,3:end);
else
    t=[];
    s=[];
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end
