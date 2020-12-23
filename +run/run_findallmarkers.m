function [M,c,cL]=run_findallmarkers(X,genelist,c)
%run R_FindAllMarkers

[c,cL]=grp2idx(c);
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/R_FindAllMarkers');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('counts.txt','file'), delete('counts.txt'); end
if exist('clusterid.txt','file'), delete('clusterid.txt'); end
if exist('genelist.txt','file'), delete('genelist.txt'); end

writematrix(X,'counts.txt');
writematrix(c,'clusterid.txt');
writematrix(genelist,'genelist.txt');

RunRcode('script.R');
if exist('output.txt','file')
    warning off
    T=readtable('output.txt');
    warning on
    n=max(T.cluster);
    M=cell(n,1);
    for k=1:n
        M{k}=string(T.gene(T.cluster==k));
    end
else
    M=[];
end
if exist('counts.txt','file'), delete('counts.txt'); end
if exist('clusterid.txt','file'), delete('clusterid.txt'); end
if exist('genelist.txt','file'), delete('genelist.txt'); end
if exist('output.txt','file'), delete('output.txt'); end  
cd(oldpth);
