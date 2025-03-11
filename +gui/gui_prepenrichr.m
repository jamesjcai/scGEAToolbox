function [outgenelist, outbackgroundlist, enrichrtype] = gui_prepenrichr(genelist, ...
    backgroundlist, questtxt, FigureHandle)

if nargin<4, FigureHandle=[]; end
enrichrtype = [];
outgenelist = [];
outbackgroundlist = [];

answer = gui.myQuestdlg(FigureHandle, questtxt);
if ~strcmp(answer, 'Yes'), return; end

answer = gui.myQuestdlg(FigureHandle, sprintf('Input list contains %d genes. Run enrichment analysis with all genes?',... 
            numel(genelist)),'',...
            {'Yes, use all genes', 'No, pick top k genes',...
            'Cancel'}, 'Yes, use all genes');
if strcmp(answer, 'Yes, use all genes')
    outgenelist = genelist;
elseif strcmp(answer, 'No, pick top k genes')
    k = gui.i_inputnumk(min([250, numel(genelist)]), 10, numel(genelist));
    if isempty(k), return; end
    outgenelist = genelist(1:k);
else
    return;
end

if ~isempty(backgroundlist)
    answer = gui.myQuestdlg(FigureHandle, 'Enrichr with background?','');
    if strcmp(answer, 'Yes')
        outbackgroundlist = backgroundlist;
    elseif strcmp(answer, 'No')
        outbackgroundlist = []; 
    else
        return;
    end
end

%enrichrtype = gui.myQuestdlg(FigureHandle, 'Select the type of Enrichr application.','', ...
%          {'Web-based', 'API-based'}, 'API-based');

enrichrtype = "Web-based";

% switch answer1 
%     case "API-based"
% 
% 
% 
%     case "Web-based"
%         fw = gui.gui_waitbar([], false, 'Sending genes to web browser...');
%         % gui.i_enrichtest(genelist, backgroundlist, numel(genelist));
%             if needbckg
%                 run.web_Enrichr_bkg(genelist, backgroundlist, numel(genelist));
%             else
%                 run.web_Enrichr(genelist, numel(genelist));
%             end
%         gui.gui_waitbar(fw, false, 'Check web browser & submit genes to Enrichr.');
%     otherwise
%         return;
% end
% 
