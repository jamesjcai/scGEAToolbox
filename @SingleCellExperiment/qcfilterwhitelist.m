function obj = qcfilterwhitelist(obj,libsize,mtratio,min_cells_nonzero,whitelist)
    if nargin<5, whitelist=[]; end
    if nargin<4 || isempty(min_cells_nonzero), min_cells_nonzero=0.01; end
    if nargin<3 || isempty(mtratio), mtratio=0.15; end
    if nargin<2 || isempty(libsize), libsize=500; end
    
    if ~isempty(whitelist)
        assert(all(ismember(whitelist,obj.g)));
        [~,idxx]=ismember(whitelist,obj.g);
        Xresv=obj.X(idxx,:);
    end
    
    [~,keptg,keptidxv]=sc_qcfilter(obj.X,obj.g,libsize,mtratio,1,...
                                   min_cells_nonzero);
    for k=1:length(keptidxv)
        obj = selectcells(obj,keptidxv{k});
    end
    [y]=ismember(obj.g,keptg);
    obj.X=obj.X(y,:);
    obj.g=obj.g(y);
    [obj.X,obj.g]=sc_rmdugenes(obj.X,obj.g);    
    if ~isempty(whitelist)
        for k=1:length(keptidxv)
            Xresv = Xresv(:,keptidxv{k});
        end
        obj.X=[obj.X; Xresv];
        obj.g=[obj.g; whitelist];
    end    

end
