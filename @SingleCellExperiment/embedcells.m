function obj = embedcells(obj,methodid,forced,usehvg)
    if nargin<4, usehvg=true; end
    if nargin<3, forced=false; end
    if nargin<2, methodid=1; end
    if isempty(obj.s) || forced
        if isstring(methodid) || ischar(methodid)
            methodid=lower(methodid);
        end
        if usehvg
            [~,X]=sc_hvg(obj.X,obj.g,true,false);
            X=X(1:min([size(X,1),2500]),:);
        else
            X=obj.X;
        end
        switch methodid
            case {1,'tsne'}
                obj.s=sc_tsne(X,3,false,true);
                obj.struct_cell_embeddings.tsne=obj.s;
            case {2,'umap'}
                obj.s=sc_umap(X,3);
                obj.struct_cell_embeddings.umap=obj.s;
            case {3,'phate'}
                obj.s=sc_phate(X,3);
                obj.struct_cell_embeddings.phate=obj.s;
        end
        disp('SCE.S added');
    else
        disp('SCE.S existed.')
        disp('Use `sce=sce.embedcells(''tSNE'',true)` to overwirte.');
    end
end
