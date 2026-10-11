function tf = i_usegpu(X, threshold)
% Return true when a CUDA GPU is available and X is large enough to benefit.
% X         - input matrix (used only for numel)
% threshold - element count above which GPU is worthwhile (default: 5e6)
%
% Every caller feeds the result to NET.PCRNET's UseGPU argument. Since
% net.pcrnet became one eigendecomposition (2026-09-24) the GPU no longer
% loses, but there is little left to win. Measured 2026-10-05 on an RTX 5060
% against 20-core BLAS, timeit, one network each:
%
%    genes x cells    CPU s    GPU s
%     1000 x   500     0.13     0.22
%     2000 x   500     0.35     0.25
%     4000 x   500     0.81     0.33
%     8000 x   500     1.93     0.82
%     2000 x  5000     1.02     0.63
%     3000 x 20000     2.57     2.28
%
% Re-measured after pcrnet's early stop (same day), which sped up both:
%
%     3000 x  2108     1.10     0.56
%     3000 x 20000     2.26     2.21
%     6000 x 20000    13.92     9.95
%
% So the GPU still pays above the threshold - 2x at the size of a
% scTenifoldXct cell type - but timeit hides its one-time cost: the first
% gpuArray call in a MATLAB session spends about 2.3 s initialising the
% device (measured inside a scTenifoldXct run). A single network in a
% fresh session is slower on the GPU; every later one is faster. The 5e6
% threshold (set 2026-08-18, when the old per-gene SVD loop ran 1-16x
% slower on the GPU) is kept.
if nargin < 2, threshold = 5e6; end
tf = gpuDeviceCount > 0 && numel(X) > threshold;
end
