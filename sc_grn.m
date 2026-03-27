function A = sc_grn(X, type, varargin)
%SC_GRN Construct single-cell gene regulatory network (scGRN)
%
%   A = sc_grn(X) uses default type 'pcrnet'.
%   A = sc_grn(X, type) where type is one of:
%       'pcrnet'          - principal component regression (parallel by default)
%       'pcrnet_batch'    - PCR (batch pagesvd)
%       'pcrnet_denoised' - tensor-denoised PCR (robust, slow)
%       'genie3'          - random forest ensemble (slow)
%       'pearson'         - thresholded Pearson correlation
%       'xicor'           - Chatterjee's xi correlation (nonlinear)
%       'distcorr'        - distance correlation (nonlinear, symmetric)
%       'mi'              - mutual information (parallel)
%       'grnformer'       - GRNFormer graph transformer (Hegde & Cheng 2026)
%
%   For 'grnformer', additional arguments are required / supported:
%       A = sc_grn(X, 'grnformer', tf_idx)
%       A = sc_grn(X, 'grnformer', tf_idx, 'GroundTruth', GT, ...)
%   where tf_idx is a logical or integer index vector of TF genes, and
%   remaining name-value pairs are forwarded to net.grnformer.
%
%   All methods are implemented in the +net/ package.
%   See also: net.pcrnet, net.genie3, net.xicornet, net.grnformer

arguments
    X {mustBeNumeric}
    type (1,1) string = "pcrnet"
end
arguments (Repeating)
    varargin
end

type = lower(string(type));

validTypes = ["pcrnet" "pcrnet_batch" "pcrnet_denoised" ...
              "genie3" "pearson" "xicor" "distcorr" "mi" "grnformer"];
if ~ismember(type, validTypes)
    error("sc_grn:InvalidType", ...
          "Type must be one of: %s", strjoin(validTypes, ", "));
end

switch type
    case "pcrnet"
        A = net.pcrnet(X);
    case "pcrnet_batch"
        A = net.pcrnet_batch(X);
    case "pcrnet_denoised"
        A = net.pcrnet_denoised(X);
    case "genie3"
        A = net.genie3(X, [], true);
    case "pearson"
        A = net.pearsonnet(X);
    case "xicor"
        A = net.xicornet(X);
    case "distcorr"
        A = net.distcorrnet(X);
    case "mi"
        A = net.minet(X);
    case "grnformer"
        if isempty(varargin)
            error("sc_grn:MissingTFIdx", ...
                "grnformer requires tf_idx as the third argument.\n" + ...
                "  Usage: sc_grn(X, 'grnformer', tf_idx, ...)");
        end
        tf_idx = varargin{1};
        extra  = varargin(2:end);
        A = net.grnformer(X, tf_idx, extra{:});
end
end
