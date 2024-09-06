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

answer = questdlg('Enrichr with background?','');
if strcmp(answer, 'Yes')
    needbckg = true;
elseif strcmp(answer, 'No')
    needbckg = false;
else
    return;
end


%answer1 = questdlg("Select Enrichr.","", ...
%     "Matlab Enrichr", "Web Enrichr", "Matlab Enrichr");
answer1 = "Web Enrichr";
switch answer1 
    case "Matlab Enrichr"



    case "Web Enrichr"
        fw = gui.gui_waitbar([], false, 'Sending genes to web browser...');
        % gui.i_enrichtest(genelist, backgroundlist, numel(genelist));
            if needbckg
                run.web_Enrichr_bkg(genelist, backgroundlist, numel(genelist));
            else
                run.web_Enrichr(genelist, numel(genelist));
            end
        gui.gui_waitbar(fw, false, 'Check web browser & submit genes to Enrichr.');
    otherwise
        return;
end

