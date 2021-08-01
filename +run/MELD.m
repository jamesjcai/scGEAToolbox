function [sample_likelihoods,T]=MELD(X,batchid)
% MELD - a graph signal processing tool used to smooth a binary variable on 
% the cell-cell graph to determine which regions of its underlying 
% data manifold are enriched or depleted in cells with a specific 
% feature.

if nargin<2
    %batchid=string([true(ceil(size(X,2)/2),1); false(floor(size(X,2)/2),1)]);
    error('USAGE: score=MELD(X,batchid)');
end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','MELD');
cd(wrkpth);

if exist('./batchid.txt','file'), delete('./batchid.txt'); end
if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output.txt','file'), delete('./output.txt'); end

X=sc_norm(X);
X=sqrt(X);

writematrix(X,'input.txt');

writematrix(batchid,'batchid.txt');

x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
disp(cmdlinestr)
[status]=system(cmdlinestr);

sample_likelihoods=[];
if status==0 && exist('output.txt','file')
    warning off
    T=readtable('output.txt',"ReadVariableNames",true);
    warning on
    sample_likelihoods=table2array(T);
end

if exist('./batchid.txt','file'), delete('./batchid.txt'); end
if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output.txt','file'), delete('./output.txt'); end
cd(oldpth);
end
