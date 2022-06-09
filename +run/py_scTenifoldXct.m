function [T]=py_scTenifoldXct(sce)

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_scTenifoldXct');
cd(wrkpth);

isdebug=true;

tmpfilelist={'X.txt','g.txt','c.txt','output.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

load(fullfile(pw1,'..','resources','Ligand_Receptor.mat'), ...
    'ligand','receptor');
validg=unique([ligand receptor]);
[y]=ismember(upper(sce.g),validg);

X=sce.X(y,:);
g=sce.g(y);

writematrix(X,'X.txt');
writematrix(g,'g.txt');
writematrix(sce.c_batch_id,'c.txt');

x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
disp(cmdlinestr)
[status]=system(cmdlinestr);

if status==0 && exist('output.txt','file')
    T=readtable('output.txt',"ReadVariableNames",false);
else
    T=[];    
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

