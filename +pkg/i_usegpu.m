function tf = i_usegpu(X, threshold)
% Return true when a CUDA GPU is available and X is large enough to benefit.
% X         - input matrix (used only for numel)
% threshold - element count above which GPU is worthwhile (default: 5e5)
if nargin < 2, threshold = 5e5; end
tf = gpuDeviceCount > 0 && numel(X) > threshold;
end
