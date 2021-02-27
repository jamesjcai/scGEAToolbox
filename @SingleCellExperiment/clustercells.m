function obj = clustercells(obj,k,methodid,forced)
    if isempty(obj.s)
        error('SCE.S is empty');
    end
    if nargin<4, forced=false; end
    if nargin<3, methodid=1; end
    if nargin<2, k=10; end
    if isempty(obj.c_cluster_id) || forced
        switch methodid
            case 1
                obj.c_cluster_id=sc_cluster_s(obj.s,k);
        end        
        disp('C_CLUSTER_ID added.');
    else
        disp('C_CLUSTER_ID existed.');
    end
end
