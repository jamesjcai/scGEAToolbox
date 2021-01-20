function obj = embedcells(obj,methodid,forced)
    if nargin<3, forced=false; end
    if nargin<2, methodid=1; end
    if isempty(obj.s) || forced
        switch methodid
            case 1
                obj.s=sc_tsne(obj.X,3,false,true);
            case 2
                obj.s=run.UMAP(obj.X,3);
            case 3
                obj.s=run.PHATE(obj.X,3);
        end        
        disp('S added.');
    else
        disp('S existed.');
    end
end
