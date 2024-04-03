function i_enrichtest(genelist, bkglist, numgene)

if nargin < 2
    bkglist = [];
end
if nargin < 3
    numgene = 200;
    k = gui.i_inputnumk(numgene, 20, 500);
    if isempty(k), return; end
    numgene = k;
end
numgene = min([length(genelist), numgene]);

run.web_Enrichr_bkg(genelist, bkglist, numgene);
