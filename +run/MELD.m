function [sample_likelihoods,T]=MELD(X,batchid,usematinput)
% MELD - a graph signal processing tool used to smooth a binary variable on 
% the cell-cell graph to determine which regions of its underlying 
% data manifold are enriched or depleted in cells with a specific 
% feature.
if nargin<3, usematinput=true; end
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
if exist('./input.mat','file'), delete('./input.mat'); end


X=sc_norm(X);
X=sqrt(X);

if usematinput
    save('input.mat','X','batchid','-v7');
else
    writematrix(X,'input.txt');
    writematrix(batchid,'batchid.txt');
end

x=pyenv;
pkg.i_add_conda_python_path;
if usematinput
    cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
else
    cmdlinestr=sprintf('"%s" "%s%sscript_csv.py"',x.Executable,wrkpth,filesep);
end
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
if exist('./input.mat','file'), delete('./input.mat'); end
cd(oldpth);
end
