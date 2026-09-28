function [modetag, usehvgs, usemarkers, markerweight] = i_validgenemode(genemode)
%I_VALIDGENEMODE Canonicalise a gene-set selection mode.
%
%   [modetag, usehvgs, usemarkers, markerweight] = ...
%       PKG.I_VALIDGENEMODE(genemode)
%
%   accepts either spelling a gene-set choice travels in and returns every
%   form of it, so all call sites agree on the names:
%
%     true  | "hvg"                    the top NUMHVG highly variable genes
%             "hvg+markers"            those plus every PanglaoDB marker in
%                                      the data, all genes equally weighted
%             "hvg+markers:weighted"   the same set, but with the markers
%                                      the HVG cut had dropped scaled up
%     false | "all"                    every gene
%
%   The logicals are the original interface of
%   SINGLECELLEXPERIMENT.EMBEDCELLS's USEHVGS argument and keep working
%   unchanged. Matching is case-insensitive.
%
%   MARKERWEIGHT is 1 for every mode but the weighted one, where it is
%   MARKERWEIGHTDEFAULT below. Only the genes the union ADDED are scaled:
%   the markers that already rank as highly variable are driving the
%   embedding anyway, and on the bundled example they are 840 of the 2000
%   HVGs, so scaling those too would rescale 42% of the feature space
%   rather than emphasise anything.
%
%   An unrecognised name raises rather than falling back to a default.
%   Reading a typo as false answers it with a 20000-gene embedding instead
%   of a 2000-gene one, and nothing about the result says the mode was not
%   understood.
%
%   See also PKG.I_GETMARKERWHITELIST, GUI.I_GETHVGNUM.

% Chosen by measurement, not taste. Scaling a gene's row by c multiplies
% its variance - and so its pull on the PCs t-SNE reduces to - by c^2. On
% the bundled 8260-cell example, mean 15-NN Jaccard against an HVG-only
% embedding in 50-PC space ran 0.788 unweighted, 0.414 at x2, 0.331 at x5
% and 0.321 at x10: past about x5 the markers have saturated the PCs and
% more weight buys nothing, so the useful band is roughly 1 to 3.
markerweightdefault = 2;

validmodes = ["hvg", "hvg+markers", "hvg+markers:weighted", "all"];

if islogical(genemode) || isnumeric(genemode)
    if ~isscalar(genemode)
        error('pkg:i_validgenemode:badMode', ...
            'A logical GENEMODE must be scalar. Got %s.', ...
            mat2str(size(genemode)));
    end
    if genemode
        modetag = "hvg";
    else
        modetag = "all";
    end
else
    modetag = lower(string(genemode));
    if ~isscalar(modetag) || ~ismember(modetag, validmodes)
        error('pkg:i_validgenemode:badMode', ...
            'GENEMODE must be true, false, or one of %s. Got "%s".', ...
            strjoin("""" + validmodes + """", ', '), ...
            strjoin(modetag, '", "'));
    end
end

usemarkers = startsWith(modetag, "hvg+markers");
usehvgs = modetag ~= "all";

markerweight = 1;
if modetag == "hvg+markers:weighted"
    markerweight = markerweightdefault;
end
end
