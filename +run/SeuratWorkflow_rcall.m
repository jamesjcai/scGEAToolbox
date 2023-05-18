function [sce]=SeuratWorkflow_rcall(X,genelist)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_SeuratWorkflow');

if ~isok, error(msg); sce=[]; return; end

if isa(X,'SingleCellExperiment') && isnumeric(genelist)
    sce=X;
    ndim=genelist;
else
    if nargin<2, error("[sce]=run.SeuratWorkflow(X,genelist)"); end
    sce=SingleCellExperiment(X,genelist);
    ndim=2;
end


if exist('./output.mat.tmp','file'), delete('./output.mat.tmp'); end
%if exist('tsneoutput.csv','file'), delete('tsneoutput.csv'); end
%if exist('umapoutput.csv','file'), delete('umapoutput.csv'); end
%if exist('activeidentoutput.csv','file'), delete('activeidentoutput.csv'); end
%sc_writefile('input.txt',sce.X,sce.g);
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./input.txt','file'), delete('./input.txt'); end
end

X=sce.X;
genelist=sce.g;
if ~iscellstr(genelist) && isstring(genelist)
    genelist=cellstr(genelist);
end
lastwarn('');
save('input.mat','X','genelist','-v6');

[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)    
    disp(warnId)
	if exist('./input.mat','file'), delete('./input.mat'); end
    disp('Writing data into input.txt...')
    sc_writefile('input.txt',sce.X,sce.g);
end

Rpath=getpref('scgeatoolbox','rexecutablepath');

if ndim==3
    pkg.RunRcode('script3d.R',Rpath);
else    
    pkg.RunRcode('script.R',Rpath);
end
%if exist('tsneoutput.csv','file')
%    T=readtable('tsneoutput.csv','ReadVariableNames',true);
%    s=table2array(T(:,2:end));
%    sce.struct_cell_embeddings.tsne=s;
%    sce.s=s;
%end
%if exist('umapoutput.csv','file')
%    T=readtable('umapoutput.csv','ReadVariableNames',true);
%    s=table2array(T(:,2:end));
%    sce.struct_cell_embeddings.umap=s;
%end
%if exist('activeidentoutput.csv','file')
%    T=readtable('activeidentoutput.csv','ReadVariableNames',true);
%    c=table2array(T(:,2:end));
%    sce.c_cluster_id=c;
%    sce.struct_cell_clusterings.seurat=c;
%    sce.c=c;
%end
if exist('output.mat','file')
    load output.mat s_tsne s_umap c_ident
    [c,~]=grp2idx(c_ident);
    sce.c_cluster_id=c;
    sce.struct_cell_clusterings.seurat=c_ident;
    sce.c=c;
    sce.struct_cell_embeddings.umap=s_umap;
    sce.struct_cell_embeddings.tsne=s_tsne;
    sce.s=s_tsne;
end
if ~isdebug
	if exist('./input.mat','file'), delete('./input.mat'); end
	if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./input.txt','file'), delete('./input.txt'); end
end
%if exist('tsneoutput.csv','file'), delete('tsneoutput.csv'); end
%if exist('umapoutput.csv','file'), delete('umapoutput.csv'); end
%if exist('activeidentoutput.csv','file'), delete('activeidentoutput.csv'); end
cd(oldpth);
end