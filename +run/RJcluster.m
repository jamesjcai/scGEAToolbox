function [s]=RJcluster(X)
%RJcluster - R pacakge
oldpth=pwd();
[isok,msg]=commoncheck_R('R_RJcluster');
if ~isok, error(msg); s=[]; return; end
if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end


[X]=sc_norm(X);
[X]=log(X+1);

writematrix(X','input.csv');

pkg.RunRcode('script.R');
if exist('output.csv','file')
    s=readmatrix('output.csv');
    s=s(:,2);
else
    s=[];
end
if exist('./input.csv','file'), delete('./input.csv'); end
if exist('./output.csv','file'), delete('./output.csv'); end
cd(oldpth);
end
