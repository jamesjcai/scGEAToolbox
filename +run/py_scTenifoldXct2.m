function [T]=py_scTenifoldXct2(sce,celltype1,celltype2,twosided)

if nargin<4, twosided=false; end

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_scTenifoldXct2');
cd(wrkpth);

isdebug=true;


idx=sce.c_cell_type_tx==celltype1 | sce.c_cell_type_tx==celltype2;
sce=sce.selectcells(idx);
sce.c_batch_id=sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx==celltype1)="Source";
sce.c_batch_id(sce.c_cell_type_tx==celltype2)="Target";
% sce=sce.qcfilter;

tmpfilelist={'X.mat','X.txt','g.txt','c.txt','output.txt', ...
             'output1.txt','output2.txt',...
             'gene_name_Source.tsv', 'gene_name_Target.tsv',...
             'pcnet_Source.npz', 'pcnet_Target.npz'};

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

% load(fullfile(pw1,'..','resources','Ligand_Receptor.mat'), ...
%     'ligand','receptor');
% validg=unique([ligand receptor]);
% [y]=ismember(upper(sce.g),validg);
% X=sce.X(y,:);
% g=sce.g(y);

% writematrix(sce.X,'X.txt');
X=sce.X;
save('X.mat','-v7.3','X');
writematrix(sce.g,'g.txt');
writematrix(sce.c_batch_id,'c.txt');

t=table(sce.g,sce.g,'VariableNames',{' ','gene_name'});
writetable(t,'gene_name_Source.tsv','filetype','text','Delimiter','\t');
writetable(t,'gene_name_Target.tsv','filetype','text','Delimiter','\t');
A1=sc_pcnetpar(sce.X(:,sce.c_cell_type_tx==celltype1));
A2=sc_pcnetpar(sce.X(:,sce.c_cell_type_tx==celltype2));
save A1.mat A1 -v7.3
save A2.mat A2 -v7.3

x=pyenv;
pkg.i_add_conda_python_path;

if ~twosided
    cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
else
    cmdlinestr=sprintf('"%s" "%s%sscript_2_sided.py"',x.Executable,wrkpth,filesep);
end
disp(cmdlinestr)
% [status]=system(cmdlinestr);
% https://www.mathworks.com/matlabcentral/answers/334076-why-does-externally-called-exe-using-the-system-command-freeze-on-the-third-call

rt=java.lang.Runtime.getRuntime(); 
% pr = rt.exec(strcat("cmd /c start /wait ", cmdlinestr));
pr = rt.exec(cmdlinestr);
[status]=pr.waitFor();

if ~twosided
    if status==0 && exist('output.txt','file')
        T=readtable('output.txt');
    end
else
    if status==0 && exist('output1.txt','file') && exist('output2.txt','file')
        T1=readtable('output1.txt');
        T2=readtable('output2.txt');
        T={T1,T2};
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
