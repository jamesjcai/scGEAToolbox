function [sce]=r_readSeuratRds(filename)

    sce=[];
    if nargin<1, error('run.r_readSeuratRds(filename)'); end
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_SeuratReadRds');
    if ~isok, error(msg); sce=[]; return; end
    isdebug=false;
    tmpfilelist={'inputrdsfile.txt','output.h5',...
        'g.csv','X.csv','umap.csv','barcodes.csv'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

    writematrix(filename,'inputrdsfile.txt');
    Rpath=getpref('scgeatoolbox','rexecutablepath');
    pkg.RunRcode('script.R',Rpath);
    
if exist('g.csv','file') && exist('output.h5','file')
    t=readtable('g.csv','Delimiter',',');
    g=string(t.x);
    X=h5read('output.h5','/X');
    sce=SingleCellExperiment(X,g);
elseif exist('g.csv','file') && exist('X.csv','file')
    t=readtable('g.csv','Delimiter',',');
    g=string(t.x);
    X=readmatrix('X.csv');    
    sce=SingleCellExperiment(X,g);    
end

if exist('barcodes.csv','file') && ~isempty(sce)
    t=readtable('barcodes.csv','Delimiter',',');
    id=string(t.x);
    sce.c_cell_id=id;
end
    
    if exist('umap.csv','file') && ~isempty(sce) && ~isempty(sce.c_cell_id)

        t=readtable('umap.csv','Delimiter',',');
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
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end