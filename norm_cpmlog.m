function [X]=norm_cpmlog(X)
    % ref: https://github.com/pachterlab/BHGP_2022/blob/main/scripts/norm_sparse.py
    % def norm_cpm_log(mtx):    mtx=sp.sparse.random(3,5,density=1)
    % https://doi.org/10.1101/2022.05.06.490859
    [X]=norm_libsize(X,1e6);
    [X]=log(X+1);
end