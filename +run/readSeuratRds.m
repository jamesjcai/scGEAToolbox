function [sce]=readSeuratRds(filename)

    sce=[];
    if nargin<1, error('run.readSeuratRds(filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratReadRds');
    if ~isok, error(msg); sce=[]; return; end
    
    if exist('./inputrdsfile.txt','file'), delete('./inputrdsfile.txt'); end
    if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./output.mat.tmp','file'), delete('./output.mat.tmp'); end
    if exist('./g.csv','file'), delete('./g.csv'); end
    if exist('./X.csv','file'), delete('./X.csv'); end    
    if exist('./umap.csv','file'), delete('./umap.csv'); end
    if exist('./barcodes.csv','file'), delete('./barcodes.csv'); end
    
    writematrix(filename,'inputrdsfile.txt');
    pkg.RunRcode('script.R');
    
if exist('g.csv','file') && exist('output.mat','file')
    t=readtable('g.csv');
    g=string(t.x);
    load output.mat X
    sce=SingleCellExperiment(X,g);
elseif exist('g.csv','file') && exist('X.csv','file')
    t=readtable('g.csv');
    g=string(t.x);
    X=readmatrix('X.csv');    
    sce=SingleCellExperiment(X,g);    
end

if exist('barcodes.csv','file') && ~isempty(sce)
    t=readtable('barcodes.csv');
    id=string(t.x);
    sce.c_cell_id=id;
end
    
%     if exist('X.csv','file')
%         t=readtable('X.csv');
%         g=string(t.Var1);
%         X=sparse(table2array(t(:,2:end)));
%         sce=SingleCellExperiment(X,g);
%     end
    if exist('umap.csv','file') && ~isempty(sce) && ~isempty(sce.c_cell_id)

        t=readtable('umap.csv');
        [y,idx]=ismember(string(t.Var1),sce.c_cell_id);
        if all(y)
            s=table2array(t(:,2:end));
            sce.s=s(idx,:);
        end
    end
%     if exist('batchid.csv','file') && ~isempty(sce)
%         t=readtable('batchid.csv');
%         id=string(t.x);
%         sce.c_batch_id=id;
%     end
%     if exist('input.txt','file'), delete('input.txt'); end
%     if exist('X.csv','file'), delete('X.csv'); end
%     if exist('umap.csv','file'), delete('umap.csv'); end
%     if exist('batchid.csv','file'), delete('batchid.csv'); end 
    if exist('./inputrdsfile.txt','file'), delete('./inputrdsfile.txt'); end
    if exist('./output.mat','file'), delete('./output.mat'); end
    if exist('./output.mat.tmp','file'), delete('./output.mat.tmp'); end    
    if exist('./g.csv','file'), delete('./g.csv'); end
    if exist('./X.csv','file'), delete('./X.csv'); end
    if exist('./umap.csv','file'), delete('./umap.csv'); end
    if exist('./barcodes.csv','file'), delete('./barcodes.csv'); end
    
    cd(oldpth);
end