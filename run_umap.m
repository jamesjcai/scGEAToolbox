function [s]=run_umap(X,plotit)
if isempty(FindRpath)
   error('Rscript.ext is not found.');
end
if nargin<2, plotit=false; end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_UMAP');
cd(pth);

if exist('output.csv','file')
    delete('output.csv');
end
%if ~exist('input.csv','file')
    csvwrite('input.csv',X');
%end
RunRcode('script.R');
if exist('output.csv','file')
    s=csvread('output.csv',1,1);
else
    s=[];
end
if exist('input.csv','file')
    delete('input.csv');
end
if exist('output.csv','file')
    delete('output.csv');
end
cd(oldpth);

if plotit && ~isempty(s)
    i_myscatter(s);
    xlabel('UMAP 1')
    ylabel('UMAP 2')
end
