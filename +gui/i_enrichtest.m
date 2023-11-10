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


% answer1 = gui.timeoutdlg(@(x){questdlg('Which functional enrichment analysis do you want to use?', 'Analysis Method', ...
%     'Enrichr', 'GOrilla', 'Enrichr+GOrilla', 'Enrichr')}, 15);
% if isempty(answer1), return; end
% switch answer1
%     case 'Enrichr'
%         run.web_Enrichr(genelist, numgene, bkglist);
%     case 'GOrilla'
%         run.web_GOrilla(genelist(1:numgene));
%     case 'Enrichr+GOrilla'
%         run.web_Enrichr(genelist, numgene, bkglist);
%         run.web_GOrilla(genelist(1:numgene));
%     otherwise
%         return;
% end
