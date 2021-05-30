function obj = clustercells(obj,k,methodid,forced)
    if isempty(obj.s)
        error('SCE.S is empty');
    end
    if nargin<4, forced=false; end
    if nargin<3, methodid='kmeans'; end
    if nargin<2, k=10; end
    if isempty(obj.c_cluster_id) || forced
        switch methodid
            case {'kmeans','snndpc'}
                id=sc_cluster_s(obj.s,k,'type',methodid);                               
            case {'sc3','simlr','soptsc','sinnlrr','specter'}
                id=sc_cluster_x(obj.X,k,'type',methodid);
        end
        obj.c_cluster_id=id(:);
        obj.struct_cell_clusterings.(methodid)=obj.c_cluster_id;
        disp('C_CLUSTER_ID added.');
    else
        disp('C_CLUSTER_ID existed.');
    end
end
