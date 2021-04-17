function obj = embedcells(obj,methodid,forced)
    if nargin<3, forced=false; end
    if nargin<2, methodid=1; end
    if isempty(obj.s) || forced
        switch methodid
            case {1,'tSNE'}
                obj.s=sc_tsne(obj.X,3,false,true);
            case {2, 'UMAP'}
                obj.s=run.UMAP(obj.X,3);
            case {3, 'PHATE'}
                obj.s=run.PHATE(obj.X,3);
        end        
        disp('SCE.S added');
    else
        disp('SCE.S existed.')
        disp('Use `sce=sce.embedcells(1,true)` to overwirte.');
    end
end
