function i_checkumapaddon()
% I_CHECKUMAPADDON - Make sure the File Exchange UMAP add-on is installed
%
%   PKG.I_CHECKUMAPADDON() returns quietly when the UMAP class from the
%   Herzenberg Lab's "Uniform Manifold Approximation and Projection (UMAP)"
%   add-on (File Exchange #71902) is on the path, and errors with install
%   instructions when it is not.
%
%   Only MATLAB releases older than R2026a need it: R2026a and newer use the
%   built-in UMAP function. The add-on used to be bundled as
%   external/ml_umap45, but File Exchange now requires a MathWorks sign-in
%   to download, so it cannot be fetched on demand; the user installs it
%   once from the Add-On Explorer instead.
%
%   See also RUN.ML_UMAP, RUN.ML_METAVIZ, PKG.I_UMAPMETHOD.

% Look the class up by exact name: EXIST('UMAP', 'file') is case-blind on
% Windows and also finds the built-in umap.m function from R2026a on.
if ~isempty(meta.class.fromName('UMAP'))
    return;
end

url = "https://www.mathworks.com/matlabcentral/fileexchange/71902";
error('scGEAToolbox:umapAddonMissing', ...
    ['UMAP on MATLAB R%s needs the "Uniform Manifold Approximation and ' ...
    'Projection (UMAP)" add-on, which is not installed. Install it from ' ...
    'the Add-On Explorer (Home > Add-Ons, search "UMAP") or from %s, ' ...
    'or upgrade to MATLAB R2026a or newer, which has UMAP built in.'], ...
    version('-release'), url);
end
