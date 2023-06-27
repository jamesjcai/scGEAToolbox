function m = min_cluster_size(label)
    m = 1000000000000;
    x = unique(label);
    for idx = 1:1:length(x) 
        value = x(idx);
        v = length(find(label==value));
        if v < m 
            m = v;
        end 
    end 
end
