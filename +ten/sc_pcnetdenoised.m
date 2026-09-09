function [A] = sc_pcnetdenoised(X, varargin)
% Construct denoised GRN using tensor decomposition on bootstrapped PCNets
%
% Thin wrapper — implementation moved to net.pcrnet_denoised.
% See also: net.pcrnet_denoised, sc_grn

[A] = net.pcrnet_denoised(X, varargin{:});
end
