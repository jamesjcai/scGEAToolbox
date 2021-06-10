function [T]=run_celltypeassignation(rankdedgenelist,k)

% A custom R code using fGSEA to assign cell type from given ranked list of
% genes
%
% see also: sc_pickmarkers
% Demo:
%gx=sc_pickmarkers(X,genelist,cluster_id,2);
%run_celltypeassignation(gx)

if nargin<2, k=5; end
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end

oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_cellTypeAssignation');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end



if exist('output.txt','file')
    delete('output.txt');
end
%if ~exist('input.txt','file')
    % txtwrite('input.txt',rankdedgenelist);
    writetable(table(rankdedgenelist),'input.txt','WriteVariableNames',false);
%end
RunRcode('script.R');
if exist('output.txt','file')
    T=readtable('output.txt');
else
    T=[];
end
T=T(1:k,:);
if exist('input.txt','file')
    delete('input.txt');
end
if exist('output.txt','file')
    delete('output.txt');
end
cd(oldpth);
end