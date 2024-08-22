function gui_enrichr(genelist, backgroundlist, questtxt)

answer = questdlg(questtxt);

if ~strcmp(answer, 'Yes'), return; end

answer = questdlg(sprintf('Input list contains %d genes. Run enrichment analysis with all genes?',... 
            numel(genelist)),'',...
            'Yes, use all genes', 'No, pick top k genes',...
            'Cancel', 'Yes, use all genes');
if strcmp(answer, 'Yes, use all genes')

elseif strcmp(answer, 'No, pick top k genes')
    k = gui.i_inputnumk(min([200, numel(genelist)]), 10, numel(genelist));
    if isempty(k), return; end
    genelist = genelist(1:k);
else
    return;
end

gui.i_enrichtest(genelist, backgroundlist, numel(genelist));

