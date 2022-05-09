function i_enrichtest(genelist,numgene)

if nargin<2, numgene=500; end
numgene=min([length(genelist), numgene]);

answer1=pkg.timeoutdlg(@(x){questdlg('Which functional enrichment analysis do you want to use?','Analysis Method', ...
    'Enrichr','GOrilla','Enrichr+GOrilla','Enrichr')},15);
if isempty(answer1), return; end
switch answer1
    case 'Enrichr'
        run.Enrichr(genelist,numgene);
    case 'GOrilla'
        run.GOrilla(genelist(1:numgene));
    case 'Enrichr+GOrilla'
        run.Enrichr(genelist,numgene);
        run.GOrilla(genelist(1:numgene));
    otherwise
        return;
end
