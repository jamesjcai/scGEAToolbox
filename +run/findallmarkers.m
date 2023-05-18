function [M,c,cL]=findallmarkers(X,genelist,c)
%run R_FindAllMarkers

[c,cL]=grp2idx(c);
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_FindAllMarkers');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

Rpath=getpref('scgeatoolbox','rexecutablepath');

[~,cmdout]=pkg.RunRcode('require.R',Rpath);
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('counts.txt','file'), delete('counts.txt'); end
if exist('clusterid.txt','file'), delete('clusterid.txt'); end
if exist('genelist.txt','file'), delete('genelist.txt'); end

if issparse(X), X=full(X); end
writematrix(X,'counts.txt');
writematrix(c,'clusterid.txt');
writematrix(genelist,'genelist.txt');

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

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
end
