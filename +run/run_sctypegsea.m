function [ct]=run_sctypegsea(X,genelist,c,species,organ)
%run scTypeGSEA
if nargin<4, species="Mouse"; end
if nargin<5, organ="Epithelium"; end
[c]=grp2idx(c);
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/R_scTypeGSEA');
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
