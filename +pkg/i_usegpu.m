function tf = i_usegpu(X, threshold)
% Return true when a CUDA GPU is available and X is large enough to benefit.
% X         - input matrix (used only for numel)
% threshold - element count above which GPU is worthwhile (default: 5e6)
%
% THRESHOLD RAISED FROM 5e5 TO 5e6 (2026-08-18). Every caller of this function
% feeds the result straight to NET.PCRNET's UseGPU argument, so the only
% question that matters is where the GPU beats the CPU for pcrnet. Measured on
% an RTX 5060 against 20-core BLAS, at 500 cells, one network each:
%
%    genes    CPU s    GPU s    GPU/CPU
%     1000      3.2     51.4     16.13x slower
%     2000     10.5    102.7      9.81x
%     4000     79.7    245.5      3.08x
%     6000    150.9    243.3      1.61x
%     8000    274.1    318.1      1.16x
%
% The GPU lost at every size tested. pcrnet's loop is n tiny svds calls on
% cells-by-n slices, so kernel launch overhead dominates, and the GPU branch
% additionally forces UseParallel=false while the CPU branch keeps
% multithreaded BLAS. The ratio does fall steadily with n and extrapolates to
% a crossover somewhere past 10,000 genes, which is why this is a raised
% threshold rather than a hard disable.
%
% The old 5e5 default meant any run past 1,000 genes at 500 cells took the
% slow path, which is exactly the regime scTenifoldNet's subsampling uses
% (TEN.I_NC calls TEN.I_XCTGRN with csubsmpl=500 and no way to opt out).
%
% NOTE the measurements vary only n, at fixed 500 cells. numel() conflates
% genes and cells, so a caller with few genes and very many cells is not
% covered by this evidence; such calls (e.g. SC_GRN on a full matrix) sit far
% above 5e6 and still take the GPU, as before.
if nargin < 2, threshold = 5e6; end
tf = gpuDeviceCount > 0 && numel(X) > threshold;
end
