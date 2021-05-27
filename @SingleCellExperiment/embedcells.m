function obj = embedcells(obj,methodtag,forced,usehvgs)
    if nargin<4 || isempty(usehvgs), usehvgs=true; end
    if nargin<3 || isempty(forced), forced=false; end
    if nargin<2, methodtag='tsne'; end
    
    % validTypes = {'tsne','umap','phate'};
    % checkType = @(x) any(validatestring(x,validTypes));
    if isempty(obj.s) || forced
        if isstring(methodtag) || ischar(methodtag)
            methodtag=lower(methodtag);
        end
        if usehvgs
            % disp('Identifying HVGs')
            [~,X]=sc_hvg(obj.X,obj.g,true,false);
            X=X(1:min([size(X,1),2000]),:);
        else
            X=obj.X;
        end
        switch methodtag
            case 'tsne'
                obj.s=sc_tsne(X,3,false,true);                
            case 'umap'
                obj.s=sc_umap(X,3);                
            case 'phate'
                obj.s=sc_phate(X,3);
        end
        obj.struct_cell_embeddings.(methodtag)=obj.s;
        disp('SCE.S added');
    else
        disp('SCE.S existed.')
        disp('Use `sce=sce.embedcells(''tSNE'',true)` to overwirte.');
    end
end
