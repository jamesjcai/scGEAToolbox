function [X,g,b]=readh5ad(fname)

isdebug=true;
if ~exist(fname,'file')
    fname
    error('FNAME is not a file.');
end

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_readh5ad');
cd(wrkpth);
tmpfilelist={'X.h5','obs.csv','var.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript.py" "%s"',x.Executable, ...
        wrkpth,filesep,fname);
disp(cmdlinestr)
[status]=system(cmdlinestr);

if status==0 && exist('obs.csv','file')
t1=readtable('obs.csv');
t2=readtable('var.csv');
g=string(t2.gene_name);
b=string(t1.index);

data=h5read('X.h5','/X/data');
indptr=h5read('X.h5','/X/indptr');
indices=h5read('X.h5','/X/indices');

%hinfo=h5info('X.h5');
%shape=double(hinfo.Groups(1).Attributes(1).Value);
%X=zeros(shape(1),shape(2));

shape = h5readatt('X.h5','/X','shape');
X=zeros(shape(1),shape(2));

% X=zeros(size(t1,1),size(t2,1));
for k=1:length(indptr)-1
    ix=indptr(k)+1:indptr(k+1);
    y=indices(ix)+1;
    X(y,k)=data(ix);
end
end
X=X.';

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
