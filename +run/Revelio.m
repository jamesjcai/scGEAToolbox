function [dc,T]=Revelio(X)
%Run Revelio
%
% see also: 
% https://github.com/danielschw188/Revelio

if isa(X,'SingleCellExperiment')
    genelist=upper(X.g);
    X=X.X;
end
if ~iscellstr(genelist) && isstring(genelist)
    genelist=cellstr(genelist);
end

oldpth=pwd();
[isok,msg]=commoncheck_R('R_Revelio');
if ~isok, error(msg); end
if exist('./input.mat','file'), delete('./input.mat'); end
if exist('./output.mat','file'), delete('./output.mat'); end
if exist('./output.csv','file'), delete('./output.csv'); end

save('input.mat','X','genelist','-v6');
pkg.RunRcode('script.R');
if exist('./output.mat','file')
    load('output.mat','dc');
end
if exist('./output.csv','file')
    T=readtable('./output.csv','ReadVariableNames',false);
end
if exist('./input.mat','file'), delete('./input.mat'); end
if exist('./output.mat','file'), delete('./output.mat'); end
if exist('./output.csv','file'), delete('./output.csv'); end
cd(oldpth);
figure;
gscatter(dc(1,:)',dc(2,:)',T.Var1)
end
