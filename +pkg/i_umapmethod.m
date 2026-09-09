function method = i_umapmethod()
% I_UMAPMETHOD - Pick a usable SGD method for external/ml_umap45
%
%   METHOD = PKG.I_UMAPMETHOD() returns 'MEX' when external/ml_umap45 ships
%   a stochastic-gradient-descent binary for the current platform, and
%   'MATLAB' when it does not. Pass the result to UMAP/setMethod.
%
%   The bundled package only has mexw64 (Windows) and mexmaci64 (Intel Mac)
%   binaries, and it does not ship the C++ sources needed to build others,
%   so Linux and Apple Silicon have no MEX to load. Its pure-MATLAB
%   implementation (optimize_layout.m) works everywhere and needs neither a
%   MEX binary nor a Java runtime, at the cost of being far slower -- on the
%   order of a hundred times on small inputs.
%
%   Without this, UMAP/setMethod('MEX') falls back to 'Java' when the binary
%   is missing, which then fails again on a MATLAB with no bundled JRE
%   (R2026b and newer) before finally landing on 'MATLAB' by way of two
%   caught exceptions and their printed reports.
%
%   Call this with external/ml_umap45 already on the path; PKG.I_ADDPATHTEMP
%   is the usual way to put it there.
%
%   See also PKG.I_ADDPATHTEMP.

persistent warned

if exist(['mexStochasticGradientDescent.' mexext], 'file')
    method = 'MEX';
    return;
end

method = 'MATLAB';
if isempty(warned)
    warned = true;
    warning('scGEAToolbox:umapNoMex', ...
        ['No UMAP MEX binary for this platform (%s); using the slower ' ...
        'pure-MATLAB implementation. Upgrade to R2026a or newer to get ' ...
        'the built-in UMAP instead.'], mexext);
end
end
