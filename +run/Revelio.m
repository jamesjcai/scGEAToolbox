function Revelio(X)
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
%if exist('./input.mat','file'), delete('./input.mat'); end
%if exist('./output.mat','file'), delete('./output.mat'); end


save('input.mat','X','genelist','-v6');
pkg.RunRcode('script.R');
%if exist('./output.mat','file')
%    load('output.mat','X','contamination')
%end
%if exist('./input.mat','file'), delete('./input.mat'); end
%if exist('./output.mat','file'), delete('./output.mat'); end
cd(oldpth);
end
