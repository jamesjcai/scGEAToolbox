function i_stabilitymethodsreport(sce, S)
%I_STABILITYMETHODSREPORT Report the resolution sweep's Methods paragraph.
%
%   pkg.i_stabilitymethodsreport(sce, S)
%
%   SC_ANNOTATIONSTABILITY writes its own Methods paragraph into
%   S.T.Properties.Description, with every number it used. This prints it
%   with its references and keeps it as SCE's "clustering" step, since the
%   sweep is how the clustering resolution was chosen. S is what
%   GUI.I_ANNOTATIONSWEEP returns.
%
%   See also SC_ANNOTATIONSTABILITY, PKG.I_METHODSREPORT.

paragraph = string(S.T.Properties.Description);
if strlength(paragraph) == 0, return; end
% The works the paragraph cites, read off the paragraph itself so the two
% cannot disagree.
cited = ["Acri et al.", "acri2026"; "Blondel et al.", "blondel2008"; ...
    "Franzén et al.", "franzen2019"; "Korsunsky et al.", "korsunsky2019"; ...
    "Reichardt and Bornholdt", "reichardt2006"];
refkeys = ["cai2020", cited(contains(paragraph, cited(:, 1)), 2)'];
pkg.i_methodsreport(sce, "clustering", paragraph, refkeys);
end
