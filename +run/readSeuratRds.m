function [sce]=readSeuratRds(filename)

    sce=[];
    if nargin<1, error('run.readSeuratRds(filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratReadRds');
    if ~isok, error(msg); end
    if exist('input.txt','file'), delete('input.txt'); end
    if exist('X.csv','file'), delete('X.csv'); end
    if exist('umap.csv','file'), delete('umap.csv'); end
    if exist('batchid.csv','file'), delete('batchid.csv'); end
    
    writematrix(filename,'input.txt');
    pkg.RunRcode('script.R');
    if exist('X.csv','file')
        t=readtable('X.csv');
        g=string(t.Var1);
        X=sparse(table2array(t(:,2:end)));
        sce=SingleCellExperiment(X,g);
    end
    if exist('umap.csv','file') && ~isempty(sce)
        t=readtable('umap.csv');
        s=table2array(t(:,2:end));
        sce.s=s;
    end    
    if exist('batchid.csv','file') && ~isempty(sce)
        t=readtable('batchid.csv');
        id=string(t.x);
        sce.c_batch_id=id;
    end
    if exist('input.txt','file'), delete('input.txt'); end
    if exist('X.csv','file'), delete('X.csv'); end
    if exist('umap.csv','file'), delete('umap.csv'); end
    if exist('batchid.csv','file'), delete('batchid.csv'); end    
    cd(oldpth);    
end