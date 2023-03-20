function obj = embedcells(obj,methodtag,forced,usehvgs, ...
                          ndim,numhvg,whitelist)
    if nargin<7, whitelist=[]; end
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
            [T,X]=sc_hvg(obj.X,obj.g,true,false,true,false,true);            
            X=X(1:min([size(X,1),numhvg]),:);
            g=T.genes(1:min([size(X,1),numhvg]));
        else
            X=obj.X;
            g=obj.g;
        end

        if ~isempty(whitelist)
            assert(all(ismember(whitelist,obj.g)));
            [~,idx]=setdiff(whitelist,g);
            if ~isempty(idx)
                [~,idxx]=ismember(whitelist,obj.g);
                Xresv=obj.X(idxx,:);
                X=[X;Xresv(idx,:)];
                g=[g;whitelist(idx)];
                %size(g)
            end
        end


        switch methodtag
            case 'tsne'
                obj.s=sc_tsne(X,ndim,false,true);                
            case 'umap'
                obj.s=sc_umap(X,ndim);                
            case 'phate'
                obj.s=sc_phate(X,ndim);
            case 'metaviz'
                obj.s=run.meta_visualization(X,ndim);
        end
        obj.struct_cell_embeddings.(methodtag)=obj.s;
        disp('SCE.S added');
    else
        disp('SCE.S existed.')
        disp('Use `sce=sce.embedcells(''tSNE'',true)` to overwirte.');
    end
end
