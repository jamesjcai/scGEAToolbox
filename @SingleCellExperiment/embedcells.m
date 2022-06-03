function obj = embedcells(obj,methodtag,forced,usehvgs,ndim,numhvg)
    if nargin<6, numhvg=2000; end
    if nargin<5, ndim=3; end
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
            [~,X]=sc_hvg(obj.X,obj.g,true,false,true,false,true);
            X=X(1:min([size(X,1),numhvg]),:);
        else
            X=obj.X;
        end

        switch methodtag
            case 'tsne'
                obj.s=sc_tsne(X,ndim,false,true);                
            case 'umap'
                obj.s=sc_umap(X,ndim);                
            case 'phate'
                obj.s=sc_phate(X,ndim);
        end
        obj.struct_cell_embeddings.(methodtag)=obj.s;
        disp('SCE.S added');
    else
        disp('SCE.S existed.')
        disp('Use `sce=sce.embedcells(''tSNE'',true)` to overwirte.');
    end
end
