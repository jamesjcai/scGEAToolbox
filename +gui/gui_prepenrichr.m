function [outgenelist, outbackgroundlist, ...
    enrichrtype] = gui_prepenrichr(genelist, ...
    backgroundlist, questtxt, parentfig)

if nargin<4, parentfig=[]; end
enrichrtype = [];
outgenelist = [];
outbackgroundlist = [];

answer = gui.myQuestdlg(parentfig, questtxt);
if ~strcmp(answer, 'Yes'), return; end

answer = gui.myQuestdlg(parentfig, ...
            sprintf(['Input list contains %d genes. ' ...
            'Run enrichment analysis with all genes?'],... 
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

    outbackgroundlist = backgroundlist;

    % answer = gui.myQuestdlg(parentfig, 'Enrichr with background?','');
    % if strcmp(answer, 'Yes')
    %     outbackgroundlist = backgroundlist;
    % elseif strcmp(answer, 'No')
    %     outbackgroundlist = []; 
    % else
    %     return;
    % end
end

%enrichrtype = gui.myQuestdlg(parentfig, 'Select the type of Enrichr application.','', ...
%          {'Web-based', 'API-based'}, 'API-based');

enrichrtype = "Web-based";

