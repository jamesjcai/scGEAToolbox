function [s]=gseapy(genelist,drdistin,usepylib)

if nargin<3, usepylib=true; end
if nargin<2, drdistin=[]; end
oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','gseapy');
cd(wrkpth);

if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end


genelist=upper(genelist);
if ~isempty(drdistin) && length(genelist)==length(drdistin)
    drdist=drdistin;    
else
    t=readtable('input_template.txt');    
    N=min([size(t,1) length(genelist)]);
    genelist=genelist(1:N);
    drdist=t.drdist(1:N);
end
T=table(genelist,drdist);
writetable(T,'input.txt','WriteVariableNames',false);
% https://gseapy.readthedocs.io/en/latest/gseapy_example.html
gsetname='KEGG_2019_Human';


    

if usepylib    
    pd = py.importlib.import_module('pandas');    
    gp = py.importlib.import_module('gseapy');
    rnk=pd.read_csv("input.txt");
    % data_mat=np.array(data_mat);
    % vars_use = py.list({py.str('KEGG_2016')});
    
    % s=np2mat(ho.Z_corr.T);
    gp.prerank(pyargs('rnk',rnk,...
        'gene_sets',py.str(gsetname),...
        'outdir',py.str(gsetname),...
        'no_plot',py.True));
    s=readtable(fullfile(gsetname,'gseapy.prerank.gene_sets.report.csv'));
    s=s(s.fdr<0.05,:);    
else    
    x=pyenv;
    cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr);
    if status==0 && exist('output.txt','file')
        s=readmatrix('output.txt');
    else
        s=[];
    end    
end
%if exist('input.txt','file'), delete('input.txt'); end
%if exist('output.txt','file'), delete('output.txt'); end
cd(oldpth);
end