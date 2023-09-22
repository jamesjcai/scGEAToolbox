function [X] = norm_cp10klog(X)
% ref: https://github.com/pachterlab/BHGP_2022/blob/main/scripts/norm_sparse.py
% def norm_cp10k_log(mtx):    mtx=sp.sparse.random(3,5,density=1)
% https://doi.org/10.1101/2022.05.06.490859
[X] = pkg.norm_libsize(X, 1e4);
[X] = log(X+1);
end