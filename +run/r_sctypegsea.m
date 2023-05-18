function [ct]=r_sctypegsea(X,genelist,c,species,organ)
%run scTypeGSEA
if nargin<4, species="Mouse"; end
if nargin<5, organ="Epithelium"; end
[c]=grp2idx(c);
if isempty(pkg.FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_scTypeGSEA');
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


pkg.RunRcode('script.R',Rpath);


if exist('output.txt','file')    
    T=readtable('output.txt','Delimiter','\t','ReadVariableNames',false);
    ct=string(T.Var1);
else
    ct=[];
end
if exist('counts.txt','file'), delete('counts.txt'); end
if exist('clusterid.txt','file'), delete('clusterid.txt'); end
if exist('genelist.txt','file'), delete('genelist.txt'); end
if exist('output.txt','file'), delete('output.txt'); end  
cd(oldpth);
end
