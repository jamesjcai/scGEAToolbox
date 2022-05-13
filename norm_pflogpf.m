function [X]=norm_pflogpf(X)

% ref: https://github.com/pachterlab/BHGP_2022/blob/main/scripts/norm_sparse.py
% https://doi.org/10.1101/2022.05.06.490859
[X]=norm_libsize(X);
[X]=log(X+1);
[X]=norm_libsize(X);

end


