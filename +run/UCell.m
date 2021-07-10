function [score]=UCell(X,genelist,tgsPos)
%UCell: Robust and scalable single-cell gene signature scoring
%UCell is an R package for calculating gene signatures
%in single-cell datasets. UCell scores, based on the 
%Mann-Whitney U statistic, are robust to dataset size 
%and heterogeneity, and their calculation demands 
%relatively less computing time and memory than other 
%available methods, enabling the processing of large 
%datasets (>10^5 cells). UCell can be applied to any
%cell vs. gene data matrix, and includes functions to 
%directly interact with Seurat objects.
% https://github.com/carmonalab/UCell

genelist=upper(genelist);
tgsPos=upper(tgsPos);
if ~iscellstr(genelist) && isstring(genelist)
    genelist=cellstr(genelist);
end
if ~iscellstr(tgsPos) && isstring(tgsPos)
    tgsPos=cellstr(tgsPos);
end
oldpth=pwd();
[isok,msg]=commoncheck_R('R_UCell');
if ~isok, error(msg); end
if exist('input.mat','file'), delete('input.mat'); end
if exist('output.csv','file'), delete('output.csv'); end

save('input.mat','X','genelist','tgsPos');
pkg.RunRcode('script.R');
if exist('output.csv','file')
    T=readtable('output.csv','ReadVariableNames',true);
    score=T.scoretype_UCell;
else
    score=[];
end
if exist('input.mat','file'), delete('input.mat'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end