function [X]=norm_pflogpf(X)
    % ref: https://github.com/pachterlab/BHGP_2022/blob/main/scripts/norm_sparse.py
    % def norm_pf_log_pf(mtx):    mtx=sp.sparse.random(3,5,density=1)
    % https://doi.org/10.1101/2022.05.06.490859
    [X]=pkg.norm_libsize(X);
    [X]=log(X+1);
    [X]=pkg.norm_libsize(X);
end